library(pheatmap)

# Read the deepTools raw read-count table.
# Run this script from the repository root.
df <- read.table(
  "data/H3K4me3_readCounts.tab",
  header = FALSE,
  sep = "\t",
  comment.char = "#"
)

# Sample names follow the column order in the input file.
colnames(df) <- c(
  "chr", "start", "end",
  "H3K4me3_Nano_Tn5_rep1",
  "H3K4me3_Nano_Tn5_rep2",
  "H3K4me3_Nano_Tn5_Beads_rep1",
  "H3K4me3_Nano_Tn5_Beads_rep2",
  "H3K4me3_yl_rep1",
  "H3K4me3_yl_rep2",
  "H3K4me3_ml_rep1",
  "H3K4me3_ml_rep2"
)

# Calculate Pearson correlations.
signal <- df[, 4:ncol(df)]
cor_matrix <- cor(signal, method = "pearson")

# Create the output directory.
dir.create("figure", showWarnings = FALSE)

# Plot and save the correlation heatmap.
pheatmap(
  cor_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  number_color = "black",
  fontsize_number = 10,
  color = colorRampPalette(c("white", "#fdbb84", "#e34a33"))(100),
  breaks = seq(0.85, 1, length.out = 101),
  border_color = NA,
  treeheight_row = 20,
  treeheight_col = 20,
  filename = "figure/H3K4me3_correlation_heatmap.pdf",
  width = 8,
  height = 8
)
