# Read input --------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
paramfile <- readLines(args[1])
query <- paramfile[1]
ref <- paramfile[2]


query_bam <- paramfile[7]
ref_bam <- paramfile[8]

outfolder <- paramfile[9]


f3_file <- paste0(outfolder, "candidates_F3_", query,"_", ref, ".tab")


library(Rsamtools)
library(GenomicAlignments)
library(dplyr)
library(GenomicRanges)
library(data.table)

# Read reads spanning the intron
f3_file <- fread(f3_file)
intron_candidates <- GRanges(unique(f3_file$coordinates_intron))


getBAMreads <- function(bamFile){
  which <- intron_candidates
  param <- ScanBamParam(scanBamFlag(isUnmappedQuery = F, 
                                  isSecondaryAlignment = F),
                      which = which)
  bam <- readGappedReads(file = bamFile, param = param, use.names = TRUE)
  GenomeInfoDb::seqlevelsStyle(bam) <-"UCSC"
  return(bam)
}

hits_search <- function(bam, ir_candidates){
  hits <- findOverlaps(query = ir_candidates, subject = GRanges(bam), type = "within")
  hits <- data.frame(hits) %>%
    group_by(queryHits) %>%
    summarise(readCount = n()) %>% 
    data.frame()
  return(hits)
}


# Get reads in BAMs

query_reads <- hits_search(getBAMreads(query_bam), intron_candidates)
ref_reads <- hits_search(getBAMreads(ref_bam), intron_candidates)


intron_candidates$query_readsCount <- 0
intron_candidates$query_readsCount[ query_reads$queryHits ] <- query_reads$readCount
intron_candidates$ref_readsCount <- 0
intron_candidates$ref_readsCount[ ref_reads$queryHits ] <- ref_reads$readCount


ir_df <- data.frame(coordinates_intron = as.character(intron_candidates), 
            query_readsCount = intron_candidates$query_readsCount,
            ref_readsCount = intron_candidates$ref_readsCount)

f3_file <- left_join(f3_file, ir_df, by = "coordinates_intron")



# Convert ids

convert_col <- function(df) {
  colnames(df) <- gsub("ref", ref, colnames(df))
  colnames(df) <- gsub("query", query, colnames(df))
  return(df)
}

f3_file <- convert_col(f3_file)


fwrite(f3_file, paste0(outfolder, "candidates_F3_", query, "_", ref, ".tab"), sep = "\t", quote = F)



