# Load necessary libraries

library(rtracklayer)
library(stringr)
library(dplyr)
library(magrittr)
library(Gviz)
library(GenomicFeatures)
library(data.table)


# Get gene symbol of each event
edb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86

get_gene_symbol <- function(coord) {
  coord_gr <- GRanges(coord)
  seqlevelsStyle(coord_gr) <- "ENSEMBL"
  gene_symb <- ensembldb::genes(edb, filter = AnnotationFilter::GRangesFilter(coord_gr))$gene_name %>%
    unique() %>%
    paste(collapse = ",")

  return(gene_symb)
}



# GENCODE Transcript annotation
txdb <- makeTxDbFromGFF(
  file = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.annotation.gff3.gz",
  organism = "Homo sapiens", dataSource = "GENCODE"
)




# Individual files ----------------------------------------------------------------------------------------------

# BAM files
drnapaths <- c(
  grep("dRNA",
    list.files("../",
      pattern = "reads_aln_sorted.bam$",
      full.names = T, recursive = T
    ),
    value = T
  ),
  list.files("../../Nanopore-WGS-data/update_Dec2020/StringTie_pipeline/isoforms_dRNA/",
    pattern = "reads_aln_sorted.bam$",
    recursive = T,
    full.names = T
  )
)
names(drnapaths) <- c(gsub(
  "\\.\\./\\/(.*)\\/dRNA.*", "\\1",
  drnapaths[-length(drnapaths)]
), "NA12878")

cdnapaths <- c(
  grep("\\/cDNA",
    list.files("../",
      pattern = "reads_aln_sorted.bam$",
      full.names = T, recursive = T
    ),
    value = T
  ),
  list.files("../../Nanopore-WGS-data/update_Dec2020/StringTie_pipeline/isoforms_cDNA/",
    pattern = "reads_aln_sorted.bam$",
    recursive = T,
    full.names = T
  )
)
names(cdnapaths) <- c(gsub(
  "\\.\\./\\/(.*)\\/cDNA.*", "\\1",
  cdnapaths[-length(cdnapaths)]
), "NA12878")

dcdnapaths <- grep("dcDNA",
  list.files("../",
    pattern = "reads_aln_sorted.bam$",
    full.names = T, recursive = T
  ),
  value = T
)
names(dcdnapaths) <- gsub(
  "\\.\\./\\/(.*)\\/dcDNA.*", "\\1",
  dcdnapaths
)


# GFF corrected files
gffcdna <- list.files("input/gff_fix/", pattern = "_cDNA.gff", full.names = T)
names(gffcdna) <- gsub(".*\\/\\/(.*)_.*\\.gff", "\\1", gffcdna)

gffdcdna <- list.files("input/gff_fix/", pattern = "_dcDNA.gff", full.names = T)
names(gffdcdna) <- gsub(".*\\/\\/(.*)_.*\\.gff", "\\1", gffdcdna)

gffdrna <- list.files("input/gff_fix/", pattern = "_dRNA.gff", full.names = T)
names(gffdrna) <- gsub(".*\\/\\/(.*)_.*\\.gff", "\\1", gffdrna)


# Main plot function ----------------------------------------------------------------------------------

candidate_plot <- function(coord_intron, coord_exon, idx, cell_line, dcDNA = F, repeat_seq = "") {
  print(idx)

  # General features

  gene_names <- get_gene_symbol(coord_intron)

  chr <- gsub("(chr.*):[0-9]+\\-[0-9]+:.", "\\1", coord_exon)
  exon_start <- as.numeric(gsub("chr.*:([0-9]+)\\-[0-9]+:.", "\\1", coord_exon))
  exon_end <- as.numeric(gsub("chr.*:[0-9]+\\-([0-9]+):.", "\\1", coord_exon))

  intron_start <- as.numeric(gsub("chr.*:([0-9]+)\\-[0-9]+:.", "\\1", coord_intron))
  intron_end <- as.numeric(gsub("chr.*:[0-9]+\\-([0-9]+):.", "\\1", coord_intron))

  # Reference

  refgenome <- GeneRegionTrack(txdb,
    chromosome = chr,
    start = exon_start,
    end = exon_end,
    name = "GENCODE v36",
    cex.title = 0.5,
    fill = "black", lwd = 0.2, col = "black",
    rotation.title = 0,
  )

  # dRNA files
  drnaReads <- AlignmentsTrack(drnapaths[cell_line],
    name = "dRNA", fill = "darkgreen",
    background.legend = "grey"
  )

  drnatxdb <- makeTxDbFromGFF(file = gffdrna[cell_line])
  drnaref <- GeneRegionTrack(drnatxdb,
    chromosome = chr,
    start = exon_start, end = exon_end,
    name = "dRNA\nTX annotation",
    col = "black",
    fill = "black",
    rotation.title = 0,
    lwd = 0.2,
    cex.title = 0.5
  )


  # cDNA files
  cdnaReads <- AlignmentsTrack(cdnapaths[cell_line],
    name = "cDNA", fill = "darkred"
  )
  cdnatxdb <- makeTxDbFromGFF(file = gffcdna[cell_line])
  cdnaref <- GeneRegionTrack(cdnatxdb,
    chromosome = chr,
    start = exon_start, end = exon_end,
    name = "cDNA\nTX annotation",
    col = "black",
    fill = "black",
    lwd = 0.2,
    cex.title = 0.5,
    rotation.title = 0
  )


  # dcDNA files( when they exist)
  if (dcDNA == T) {
    dcdnaReads <- AlignmentsTrack(dcdnapaths[cell_line],
      name = "dcDNA", fill = "darkgoldenrod",
    )

    dcdnatxdb <- makeTxDbFromGFF(file = gffdcdna[cell_line])
    dcdnaref <- GeneRegionTrack(dcdnatxdb,
      chromosome = chr,
      start = exon_start, end = exon_end,
      name = "dcDNA\nTX annotation",
      fill = "black", lwd = 0.2,
      col = "black", fill = "black",
      cex.title = 0.5
    )

    ht <- HighlightTrack(
      trackList = list(refgenome, cdnaReads, cdnaref, dcdnaReads, dcdnaref, drnaReads, drnaref),
      start = intron_start,
      end = intron_end,
      chromosome = chr,
      alpha = 0.8,
      lwd = 0.01
    )
  }
  else {
    ht <- HighlightTrack(
      trackList = list(refgenome, cdnaReads, cdnaref, drnaReads, drnaref),
      start = intron_start,
      end = intron_end,
      chromosome = chr,
      alpha = 0.8,
      lwd = 0.01
    )
  }


  # Gviz Plot

  if (repeat_seq != "") {
    main_title <- paste0(coord_intron, " ", gene_names, "\nCell line: ", cell_line, "\n", "Direct repeat sequence: ", repeat_seq)
  }
  else {
    main_title <- paste0(coord_intron, " ", gene_names, "\nCell line: ", cell_line)
  }
  plotTracks(ht,
    title.width = 1,
    cex.main = 1.5,
    chromosome = chr,
    from = exon_start,
    to = exon_end,
    type = c("coverage", "sashimi", "pileup"),
    sashimiScore = 5,
    lwd.sashimiMax = 1,
    sashimiNumbers = TRUE,
    extend.right = 200, extend.left = 200,
    min.height = 0,
    coverageHeight = 0.08,
    sashimiHeight = 0.05,
    minCoverageHeight = 0,
    lwd.coverage = NA,
    transcriptAnnotation = "gene",
    fill.reads = "darkblue",
    lwd.reads = 0.01,
    lwd.gap = 0.2,
    background.title = "white",
    col.frame = "white",
    col.axis = "grey80",
    col = "grey80",
    showAxis = F,
    col.title = "black",
    col.sashimi = "black",
    fill.sashimi = "white",
    fontcolor.legend = "black",
    cex.sashimi = 0.6,
    background.legend = "black",
    main = main_title,
    # cex.title = 0.5
    rotation.title = 0
  )
}

# Test -----------------------------------------------------------------------------------------------
candidate_plot("chr16:28932452-28932581:+", "chr16:28932346-28932612:+", 1, "NA12878", dcDNA = F)



# Plot all the intron candidates ---------------------------------------------------------------------
library(tidyr)
jx_ratios_red <- jx_ratios %>%
  pivot_wider(names_from = lib_dt, values_from = c(readsCount, jx, jx_ratio)) %>%
  mutate(cdna_dcdna = ifelse(!is.na(readsCount_dcDNA), T, F)) %>%
  mutate(repeat_seq = ifelse(kmer_len == 3, "", repeat_motif))


pdf("pdf_plots/SuppFile_candidates_identified.pdf", width = 10, height = 15, useDingbats = F, family = "ArialMT")
mapply(
  candidate_plot, jx_ratios_red$coordinates_intron, jx_ratios_red$coordinates_exon,
  seq(1, nrow(jx_ratios_red)), jx_ratios_red$cell_line,
  jx_ratios_red$cdna_dcdna,
  jx_ratios_red$repeat_seq
)
dev.off()
