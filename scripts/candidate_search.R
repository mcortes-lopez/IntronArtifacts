library(data.table)
library(tidyr)
library(dplyr)

readformat <- function(file, fclass) {
  fclass_o <- ifelse(fclass == "ref", "query", "ref")

  tracking_file <- fread(file, header = F)
  colnames(tracking_file) <- c("transfrag_id", "locus_id", "reference_gene_id", "class_code", "extra")

  tracking_file <- tracking_file %>%
    separate(extra, paste0(c("qJ", "gene_id", "transcript_id", "num_exons", "FPKM", "TPM", "cov", "len"), "_", fclass), sep = ":|\\|") %>%
    separate(reference_gene_id, paste0(c("gene_id", "transcript_id"), "_", fclass_o), sep = "\\|")

  tonum <- paste0(c("num_exons", "FPKM", "TPM", "cov", "len"), "_", fclass)
  tracking_file <- tracking_file %>%
    mutate_at(vars({{ tonum }}), as.numeric)
}


# Read input --------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
paramfile <- readLines(args[1])
query <- paramfile[1]
ref <- paramfile[2]

query_ref <- readformat(paramfile[3], "ref")
ref_query <- readformat(paramfile[4], "query")

head(query_ref)
head(ref_query)

query_gff <- paramfile[5]
ref_gff <- paramfile[6]

query_bam <- paramfile[7]
ref_bam <- paramfile[8]

outfolder <- paramfile[9]



write_convert_col <- function(df, file_id) {
  colnames(df) <- gsub("ref", ref, colnames(df))
  colnames(df) <- gsub("query", query, colnames(df))
  fwrite(df, paste0(outfolder, file_id, query, "_", ref, ".tab"), sep = "\t", quote = F)
}




# Import gffcompare data ------------------------------------------------------------------------------------------------------------

# Get the TPM information (ref)
ref_tx_info <- query_ref %>%
  select("gene_id_ref", "transcript_id_ref", "FPKM_ref", "TPM_ref", "cov_ref", "len_ref") %>%
  unique()


ref_query.diff <- ref_query %>%
  filter(class_code != "=")

colnames(ref_query.diff)
colnames(ref_tx_info)

# Join with the previously obtained ref info
rec_candidates <- merge(ref_tx_info, ref_query.diff, by = c("gene_id_ref", "transcript_id_ref"))

# Save candidates from the first search
write_convert_col(rec_candidates, "candidates_F0_")

# Intron in exon search -----------------------------------------------------------------------------------

# The gffs have to be parsed first with the following line:
# grep -v -P '\t\.\t\.\t' original.gff > data/new.gff

file_pfx <- gsub("output\\/", "", outfolder)

system(paste0("grep -v -P \'", "\t", "\\", ".", "\t", "\\", ".", "\t\' ", ref_gff, " > ", file_pfx, "_ref_tmp.gff"))
system(paste0("grep -v -P \'", "\t", "\\", ".", "\t", "\\", ".", "\t\' ", query_gff, " >", file_pfx, "_query_tmp.gff"))


library(GenomicFeatures)


# Create transcript DB from gff
query.tx <- makeTxDbFromGFF(paste0(file_pfx, "_query_tmp.gff"))
ref.tx <- makeTxDbFromGFF(paste0(file_pfx, "_ref_tmp.gff"))


ref_exons <- exonsBy(ref.tx, by = "tx", use.names = T) # exons from the directRNA set
query_introns <- intronsByTranscript(query.tx, use.names = T) # introns from query, saved in a list with the TX name



# Search in the J events the "introns" that are fully contained in the exons of "ref"

contained_introns <- function(ref_txid, query_txid) {
  if (query_txid %in% names(query_introns) & ref_txid %in% names(ref_exons)) {
    ovp <- findOverlapPairs(query_introns[[query_txid]], ref_exons[[ref_txid]], type = "within")
    if (length(ovp) > 0) { # If at least one intron is contained
      coords <- as.character(second(ovp)) # Get coordinates of the ref exon
      coords_intron <- as.character(first(ovp)) # Get coordinates of the query intron
      return(list(
        "n_candidates" = length(ovp),
        "coordinates_exon" = coords,
        "coordinates_intron" = coords_intron
      )) # Return N of introns and exon location
    }
  }
}


# From the J events apply the search function
candidate_ex <- mapply(contained_introns,
  rec_candidates$transcript_id_ref,
  rec_candidates$transcript_id_query,
  USE.NAMES = F
)

# Name the list with the ids of the TXs
names(candidate_ex) <- mapply(function(x, y) {
  paste0(x, "_", y)
}, rec_candidates$transcript_id_ref, rec_candidates$transcript_id_query)

# Transform to DF, NULL records eliminated
candidate_ex <- candidate_ex %>%
  rbindlist(use.names = T, idcol = "ref_query", fill = T) %>%
  separate(ref_query, c("transcript_id_ref", "transcript_id_query"), sep = "_")


# Join to get the info
candidate_ex.tab <- merge(candidate_ex, rec_candidates,
  by = c("transcript_id_ref", "transcript_id_query")
)


# Save candidates from the first search
write_convert_col(candidate_ex.tab, "candidates_F1_1_")



## Search for exons type 2 -----------------------------------------------------------------------------------------------------

ref_introns <- intronsByTranscript(ref.tx, use.names = T)



### Candidate retained introns between two exons

contained_introns_two <- function(ref_txid, query_txid, idx) {
  print(idx)
  if (query_txid %in% names(query_introns) & ref_txid %in% names(ref_exons)) {
    if (length(ref_exons[[ref_txid]]) >= 2) {
      # Only search in transcripts starting from 2 exons
      int_filt <- subsetByOverlaps(query_introns[[query_txid]], ref_introns[[ref_txid]], type = "within", invert = T)

      event_exons <- ref_exons[[ref_txid]]

      # Get info on introns that overlap with exons
      ov_events <- mergeByOverlaps(event_exons, int_filt, type = "any")


      if (nrow(ov_events) > 1) { # 1 intron should overlap at least 2 exons
        if (any(table(ov_events$int_filt) == 2)) {
          # Count the overlaps
          int_counts <- as.data.frame(table(ov_events$int_filt))
          int_counts <- subset(int_counts, Freq == 2) # filter by 2 exons overlap

          ov_events <- subset(ov_events, as.character(ov_events$int_filt) %in% int_counts$Var1)

          # Separate introns to evaluate each individually
          sep_df <- split.data.frame(ov_events, ov_events$int_filt)

          exon_intron_list <- lapply(sep_df, function(x) {
            int_ev <- unique(x$int_filt)
            rank_diff <- x$exon_rank[2] - x$exon_rank[1] # Overlapping exons are continuos
            # print(x$exon_rank)
            upstream_e_st <- start(x$event_exons[1])
            downstream_e_end <- end(x$event_exons[2])


            # The exons should flank the intron, not being contained on it
            if (rank_diff == 1 &
              start(int_ev) > upstream_e_st &
              end(int_ev) < downstream_e_end) {

              # Exons start-end
              exon_coords <- x$event_exons[1]
              end(exon_coords) <- end(x$event_exons)[2]
              exon_coords <- as.character(exon_coords)


              intron_coords <- as.character(int_ev)

              # print(exon_coords)

              return(data.frame("coordinates_exon" = exon_coords, "coordinates_intron" = intron_coords))
            }
          })
        }
      }
    }
  }
}


candidate_ex2 <- mapply(contained_introns_two,
  rec_candidates$transcript_id_ref,
  rec_candidates$transcript_id_query,
  1:nrow(rec_candidates),
  USE.NAMES = F
)

names(candidate_ex2) <- mapply(
  function(x, y) {
    paste0(x, "_", y)
  },
  rec_candidates$transcript_id_ref,
  rec_candidates$transcript_id_query
)


candidate_ex_type2 <- lapply(
  candidate_ex2[sapply(candidate_ex2, length) == 1],
  function(x) rbindlist(x, fill = T)
) %>%
  rbindlist(use.names = T, idcol = "ref_query", fill = T) %>%
  separate(ref_query, c("transcript_id_ref", "transcript_id_query"), sep = "_")


candidate_ex_type2.tab <- merge(candidate_ex_type2,
  rec_candidates,
  by = c("transcript_id_ref", "transcript_id_query")
)

# Save candidates Filter 1.2
write_convert_col(candidate_ex_type2.tab, "candidates_F1_2_")


# Merge F 1.1 and 1.2
candidates_f1 <- rbindlist(list("type_1" = candidate_ex.tab, "type_2" = candidate_ex_type2.tab), fill = T, idcol = "type")

# Filter by intronic reads absent in the direct RNA seq data set ---------------------------------------------------
library(Rsamtools)
library(GenomicAlignments)

get_junctions <- function(coords, bamfile) {
  which <- GRanges(coords)
  param <- ScanBamParam(which = which)
  aln <- readGAlignments(bamfile, use.names = TRUE, param = param)
  jx <- summarizeJunctions(aln)
  return(jx)
}

intron_filter <- function(intron_coords, exon_coords) {
  jx_query <- get_junctions(exon_coords, query_bam)
  jx_ref <- get_junctions(exon_coords, ref_bam)
  intron_q <- GRanges(intron_coords)

  jx_query <- subsetByOverlaps(jx_query, intron_q, type = "equal", ignore.strand = T)$score
  jx_ref <- subsetByOverlaps(jx_ref, intron_q, type = "equal", ignore.strand = T)$score

  return(list("query_jx" = jx_query, "ref_jx" = jx_ref))
}



jx_quant <- mapply(intron_filter, candidates_f1$coordinates_intron, candidates_f1$coordinates_exon, SIMPLIFY = F)

jx_count <- data.table::rbindlist(jx_quant, use.names = T, fill = T)


# Add junction counts
candidates_f1 <- cbind(candidates_f1, jx_count)

# Transform NA
candidates_f1 <- candidates_f1 %>%
  mutate_at(vars(contains("jx")), ~ replace(., is.na(.), 0))

candidates_f2 <- candidates_f1 %>%
  filter(!ref_jx > 0) %>% # Filter all the junctions with reads in ref
  filter(query_jx > 4) %>% # Filter for at least 5 reads in the junctions
  filter(cov_query >= 5 & cov_ref >= 5) # Filter by coverage

candidates_f2 <- candidates_f2[!duplicated(candidates_f2$coordinates_intron), ]

# Save candidates
write_convert_col(candidates_f2, "candidates_F2_")


# Splice site search -------------------------------------------------------------------------------------------------

intron_candidates <- candidates_f2

intron.gr <- GRanges(intron_candidates$coordinates_intron)

ss_r <- resize(intron.gr, width = 2, fix = "start")
ss_l <- resize(intron.gr, width = 2, fix = "end")

library(BSgenome.Hsapiens.NCBI.GRCh38)
seqlevelsStyle(Hsapiens) <- "UCSC"

ss_df <- data.frame(
  "ss_r" = getSeq(Hsapiens, ss_r, as.character = T),
  "ss_l" = getSeq(Hsapiens, ss_l, as.character = T)
)


intron_candidates <- cbind(intron_candidates, ss_df)

intron_candidates <- intron_candidates %>%
  filter(!ref_jx > 0) %>% # Filter all the junctions with reads in ref
  filter(query_jx > 4) %>% # Filter for at least 5 reads in the junctions
  filter(cov_query >= 5 & cov_ref >= 5) %>% # Filter by coverage
  filter(!(ss_r == "GT" & ss_l == "AG"))


non_ss_intron_candidates <- intron_candidates[!duplicated(intron_candidates$coordinates_intron), ]

# Save candidates
write_convert_col(non_ss_intron_candidates, "candidates_F3_")

# Get sequences
seq2search_100nt <- getSeq(Hsapiens, GRanges(non_ss_intron_candidates$coordinates_intron) + 50)
names(seq2search_100nt) <- non_ss_intron_candidates$coordinates_intron

# writeLines(seq2search_100nt, "sequences_to_predict.txt")
# Sequences saved as:
# 50+intron+50


# Save candidates
ShortRead::writeFasta(
  seq2search_100nt,
  paste0(outfolder, "candidates_F3_seq_p50_", query, "_", ref, ".tab")
)


system(paste0("rm ", file_pfx, "_query_tmp.gff ", file_pfx, "_ref_tmp.gff"))


# Reproducibility
# writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
