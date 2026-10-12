
# ============================================================
# Peak Overlap Analysis in Rice
# Comparison of Nano-Tn5, Nano-Tn5-Beads, and Bulk ChIP
# ============================================================

# Step 1. Load required libraries
library(GenomicRanges)
library(ChIPpeakAnno)

# Step 2. Set input and output directories
data_dir <- "MACS"
output_dir <- "results"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Step 3. Define peak files
peak_files <- list(
  Nano_Tn5_rep1 = "5/NA_peaks.broadPeak",
  Nano_Tn5_rep2 = "6/NA_peaks.broadPeak",
  Nano_Tn5_Beads_rep1 = "7/NA_peaks.broadPeak",
  Nano_Tn5_Beads_rep2 = "8/NA_peaks.broadPeak",
  Bulk_ChIP_yl_rep1 = "1611.broadPeak",
  Bulk_ChIP_yl_rep2 = "1612.broadPeak",
  Bulk_ChIP_ml_rep1 = "1645.broadPeak",
  Bulk_ChIP_ml_rep2 = "1646.broadPeak"
)

# Step 4. Read broadPeak files and convert to GRanges
read_broadpeak <- function(file) {
  peak <- read.table(
    file.path(data_dir, file),
    header = FALSE,
    col.names = c(
      "chrom", "start", "end", "name", "score",
      "strand", "signal", "pval", "qval"
    )
  )

  # Convert BED 0-based starts to 1-based coordinates
  peak$start <- peak$start + 1

  makeGRangesFromDataFrame(
    peak,
    seqnames.field = "chrom",
    start.field = "start",
    end.field = "end",
    ignore.strand = TRUE
  )
}

peak_list <- lapply(peak_files, read_broadpeak)

# Step 5. Merge peaks from biological replicates
peaks_group1 <- reduce(c(
  peak_list$Nano_Tn5_rep1,
  peak_list$Nano_Tn5_rep2
))

peaks_group2 <- reduce(c(
  peak_list$Nano_Tn5_Beads_rep1,
  peak_list$Nano_Tn5_Beads_rep2
))

peaks_group3 <- reduce(c(
  peak_list$Bulk_ChIP_yl_rep1,
  peak_list$Bulk_ChIP_yl_rep2
))

peaks_group4 <- reduce(c(
  peak_list$Bulk_ChIP_ml_rep1,
  peak_list$Bulk_ChIP_ml_rep2
))

# Step 6. Identify overlapping peaks
overlap_result <- findOverlapsOfPeaks(
  group1 = peaks_group1,
  group2 = peaks_group2,
  group3 = peaks_group3,
  group4 = peaks_group4
)

# Step 7. Define colors and sample names
colors <- c(
  "#264653",
  "#2A9D8F",
  "#E9C46A",
  "#F4A261"
)

sample_names <- c(
  "snCT_Nano-Tn5",
  "snCT_Nano-Tn5_Beads",
  "Bulk_ChIP_yl",
  "Bulk_ChIP_ml"
)

# Step 8. Plot and save Venn diagram
pdf(
  file.path(output_dir, "venn_diagram_peaks_rice.pdf"),
  width = 15,
  height = 15
)

makeVennDiagram(
  overlap_result,
  NameOfPeaks = sample_names,
  fill = colors,
  lty = "blank",
  main = "Overlap of Peaks",
  cat.col = colors,
  cat.cex = 1.2,
  cex = 1.2,
  main.cex = 1.5
)

dev.off()
