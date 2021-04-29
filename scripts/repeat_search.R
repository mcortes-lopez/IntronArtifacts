## Analysis of repeat presence in candidates

library(seqinr)
library(Biostrings)
library(data.table)
library(dplyr)
library(tidyr)


# Get the sequence around the splice sites


# Splice site search ----------------------------------------------------------------
library(BSgenome.Hsapiens.NCBI.GRCh38)
Hsapiens
seqlevelsStyle(Hsapiens) <- "UCSC"
library(GenomicRanges)
intron.gr <- unique(GRanges(fread("output/candidate_introns.tab")$coordinates_intron))
names(intron.gr) <- as.character(intron.gr)

ss_r <- resize(intron.gr, width = 2, fix = "start")
ss2_r <- getSeq(Hsapiens, ss_r) # 5 SS Dinucleotide sequence
ss_r <- resize(ss_r, width = 20, fix = "center")

ss_l <- resize(intron.gr, width = 2, fix = "end")
ss2_l <- getSeq(Hsapiens, ss_l) # 3 SS Dinucleotide sequence
ss_l <- resize(ss_l, width = 20, fix = "center")


ss_r_seq <- getSeq(Hsapiens, ss_r)
ss_l_seq <- getSeq(Hsapiens, ss_l)

# names(ss_r_seq) <- names(seqp50)
# names(ss_l_seq) <- names(seqp50)

# ShortRead::writeFasta(ss_r_seq, "output/candidate_splice_sites_5ss.fa")
# ShortRead::writeFasta(ss_l_seq, "output/candidate_splice_sites_3ss.fa")


# Search matching repeats --------------------------------------------

# Function to search for the splice site overlap
match_ss <- function(st_rep, end_rep, sstype) {
  pos_in_rep <- seq(st_rep, end_rep)
  if (sstype == "3ss") {
    ss_match <- any(pos_in_rep %in% seq(8, 10))
  }
  else {
    ss_match <- any(pos_in_rep %in% seq(10, 12))
  }
  return(ss_match)
}

# Main repeat search function
repeat_search <- function(ss5, ss3, ss5_seq, ss3_seq, mlim_km = 3, lim_km = 14, idx) {
  print(idx)
  ss5 <- as.character(ss5)
  ss3 <- as.character(ss3)


  kmer_search <- lapply(seq(mlim_km, lim_km), function(klen) {
    # Separate the 20nt into different lenght kmers
    kmers <- substring(ss5, first = 1:(nchar(ss5) - (klen - 1)), last = klen:nchar(ss5))

    names(kmers) <- kmers
    kmers_seq <- DNAStringSet(x = kmers, use.names = T)

    matchkmers_ss5 <- matchPDict(pdict = kmers_seq, DNAString(ss5), max.mismatch = 0, with.indels = F)
    matchkmers_ss3 <- matchPDict(kmers_seq, DNAString(ss3), max.mismatch = 0, with.indels = F) # Search in the other splice site
    matchkmers_ss5 <- unlist(matchkmers_ss5)
    matchkmers_ss3 <- unlist(matchkmers_ss3)
    matchkmers_ss3 <- as.data.frame(matchkmers_ss3)
    matchkmers_ss5 <- as.data.frame(matchkmers_ss5)

    # Data frame with repeat info
    matchkmers <- merge(matchkmers_ss5, matchkmers_ss3,
      by = c("names", "width"),
      suffixes = c(".ss5", ".ss3")
    )

    # Only return kmers with a match in both splice sites
    if (nrow(matchkmers) >= 1) {
      matchkmers <- matchkmers %>%
        rowwise() %>%
        mutate(
          ss_match_5 = match_ss(start.ss5, end.ss5, "5ss"),
          ss_match_3 = match_ss(start.ss3, end.ss3, "3ss")
        ) %>%
        ungroup() %>%
        filter(ss_match_5 == T & ss_match_3 == T)

      return(matchkmers)
    }
  })

  names(kmer_search) <- seq(mlim_km, lim_km)
  kmer_res <- rbindlist(kmer_search, idcol = "kmer_len")

  kmer_res$ss5_seq <- as.character(ss5_seq)
  kmer_res$ss3_seq <- as.character(ss3_seq)

  return(kmer_res)
}

# For multiple core processing
repeat_ss <- parallel::mcmapply(repeat_search, ss_l_seq, ss_r_seq, ss2_l, ss2_r, 3, 20,
  seq_along(ss_l_seq),
  SIMPLIFY = F, mc.cores = 12
)

# Single core
# repeat_ss <- mapply(repeat_search, ss_l_seq, ss_r_seq, ss2_l, ss2_r, 3,20,
#                                seq_along(ss_l_seq), SIMPLIFY = F)

names(repeat_ss) <- as.character(intron.gr)

repeat_ss <- rbindlist(repeat_ss, idcol = "event_coord", use.names = T, fill = T)

repeat_ss$kmer_len <- as.numeric(repeat_ss$kmer_len)


longest_repeats <- repeat_ss %>%
  group_by(event_coord) %>%
  arrange(desc(kmer_len)) %>%
  dplyr::slice(1) %>%
  ungroup()

all_introns <- data.frame(
  event_coord = as.character(intron.gr),
  ss5_seq = as.character(ss2_l),
  ss3_seq = as.character(ss2_r)
)

# Complement the information with the introns that do not have repeats

repeats_match_ss <- left_join(all_introns, longest_repeats, by = c("event_coord", "ss5_seq", "ss3_seq"))

repeats_match_ss[is.na(repeats_match_ss$kmer_len), ]$kmer_len <- 0

# fwrite(repeats_match_ss, "output/candidates_longest_repeats.tab")
# fwrite(repeat_ss, "output/candidates_direct_repeats_full.tab")
