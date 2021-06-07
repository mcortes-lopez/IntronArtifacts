library(data.table)
library(dplyr)
library(rtracklayer)
library(GenomicFeatures)
library(GenomicRanges)


event_candidates <- unique(GenomicRanges::GRanges(readLines("bash/intron_candidates.list")))


# Get the overlapping events ------------------------------------------------------------------
find_event_overlap <- function(dataset_gr, min_diff = 5){
 
 ovp <- findOverlapPairs(event_candidates, dataset_gr)
  
  st_diff <- S4Vectors::start(S4Vectors::first(ovp)) - S4Vectors::start(S4Vectors::second(ovp))
  end_diff <- S4Vectors::end(S4Vectors::first(ovp)) - S4Vectors::end(S4Vectors::second(ovp))
  
  potential_ov <- ovp[abs(end_diff)<= min_diff & abs(st_diff) <= min_diff]

  ov_events <- unique(as.character(S4Vectors::first(potential_ov)))
  
  return(ov_events)
}

# Sequel II UHRR dataset 
sequel_ii <- import.bed("input/data_external/introns_Sequel_II_UHRR.bed")
writeLines(find_event_overlap(sequel_ii),  "output/artifacts_found_in_PACBIO_Sequel_II_UHRR.list")

# Alzheimer dataset 
alzheimer_isoseq_pre <- fread("input/data_external/Alzheimer_IsoSeq2019.preFilter.sqanti_junctions.txt")
colnames(alzheimer_isoseq_pre) <- gsub("genomic_|_coord","", colnames(alzheimer_isoseq_pre))
alzheimer_isoseq_pre.gr <-GRanges(alzheimer_isoseq_pre)
writeLines(find_event_overlap(alzheimer_isoseq_pre.gr),  "output/artifacts_found_in_PACBIO_Sequel_II_Alzheimer.prefilter.list")


alzheimer_isoseq_pos <- fread("input/data_external/Alzheimer_IsoSeq2019.postFilter.sqanti_junctions.txt")
colnames(alzheimer_isoseq_pos) <- gsub("genomic_|_coord","", colnames(alzheimer_isoseq_pos))
alzheimer_isoseq_pos.gr <-GRanges(alzheimer_isoseq_pos)
writeLines(find_event_overlap(alzheimer_isoseq_pos.gr),  "output/artifacts_found_in_PACBIO_Sequel_II_Alzheimer.postfilter.list")


# Melanoma dataset 
melanoma_isoseq_pre <- fread("input/data_external/final.rep_junctions.txt")
colnames(melanoma_isoseq_pre) <- gsub("genomic_|_coord","", colnames(melanoma_isoseq_pre))
melanoma_isoseq_pre.gr <-GRanges(melanoma_isoseq_pre)
writeLines(find_event_overlap(melanoma_isoseq_pre.gr),  "output/artifacts_found_in_PACBIO_Sequel_II_Melanoma.prefilter.list")


melanoma_isoseq_pos <- fread("input/data_external/final_filtered_junctions.txt")
colnames(melanoma_isoseq_pos) <- gsub("genomic_|_coord","", colnames(melanoma_isoseq_pos))
melanoma_isoseq_pos.gr <-GRanges(melanoma_isoseq_pos)
writeLines(find_event_overlap(melanoma_isoseq_pos.gr),  "output/artifacts_found_in_PACBIO_Sequel_II_Melanoma.postfilter.list")




