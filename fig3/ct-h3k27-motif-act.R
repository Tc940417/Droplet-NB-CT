
# ============================================================
# chromVAR Motif Activity Analysis in Rice
# H3K27ac Motif Deviation Heatmap and Top Motif Visualization
# ============================================================

# Step 1. Load required libraries
library(Seurat)
library(Signac)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)

# Step 2. Set input and output directories
data_dir <- "data"
output_dir <- "results/chromVAR"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Step 3. Read Seurat object
h3k27 <- readRDS(
  file.path(data_dir, "ct.h3k27.rds")
)

# Step 4. Extract chromVAR deviation scores
motif_matrix <- GetAssayData(
  h3k27,
  assay = "chromvar",
  layer = "data"
)

# Step 5. Calculate mean motif deviation per cell type
cell_types <- h3k27$cell_type
names(cell_types) <- colnames(h3k27)

# Match cells and exclude missing annotations
valid_cells <- colnames(motif_matrix)[
  !is.na(cell_types[colnames(motif_matrix)])
]

motif_matrix <- motif_matrix[, valid_cells, drop = FALSE]
cell_types <- cell_types[valid_cells]

motif_avg_matrix <- sapply(
  unique(as.character(cell_types)),
  function(ct) {
    cells <- names(cell_types)[as.character(cell_types) == ct]
    Matrix::rowMeans(motif_matrix[, cells, drop = FALSE])
  }
)

motif_avg_matrix <- as.matrix(motif_avg_matrix)

# Step 6. Plot chromVAR motif deviation heatmap
ht <- Heatmap(
  motif_avg_matrix,
  name = "Deviation",
  col = colorRamp2(
    c(-2, 0, 2),
    c("#000A3D", "white", "#D7301F")
  ),
  show_row_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_title = "H3K27ac Motif Activity in Rice"
)

# Step 7. Save heatmap
pdf(
  file.path(output_dir, "H3K27ac_chromVAR_heatmap.pdf"),
  width = 8,
  height = 10
)

draw(ht)
dev.off()

# Step 8. Identify top 10 motifs for each cell type
top10_list <- lapply(
  colnames(motif_avg_matrix),
  function(ct) {

    vals <- motif_avg_matrix[, ct]
    vals <- sort(vals[is.finite(vals)], decreasing = TRUE)

    top_motifs <- head(names(vals), 10)

    data.frame(
      cell_type = ct,
      motif = top_motifs,
      deviation = vals[top_motifs],
      row.names = NULL
    )
  }
)

top10_motif_df <- do.call(rbind, top10_list)

# Step 9. Save average deviation matrix and top motifs
write.csv(
  motif_avg_matrix,
  file = file.path(output_dir, "H3K27ac_motif_avg_matrix.csv"),
  row.names = TRUE,
  quote = FALSE
)

write.csv(
  top10_motif_df,
  file = file.path(output_dir, "H3K27ac_top10_motifs.csv"),
  row.names = FALSE,
  quote = FALSE
)

# Step 10. Plot motif sequence logos
top_motifs <- unique(top10_motif_df$motif)

# Retain motifs available in the ATAC motif object
available_motifs <- colnames(Motifs(h3k27[["ATAC"]]))
top_motifs <- intersect(top_motifs, available_motifs)

if (length(top_motifs) > 0) {

  p <- MotifPlot(
    object = h3k27,
    motifs = top_motifs,
    assay = "ATAC"
  )

  ggsave(
    filename = file.path(output_dir, "H3K27ac_top_motif_logos.pdf"),
    plot = p,
    device = "pdf",
    width = 12,
    height = 10
  )

  print(p)

} else {
  warning("No matching motifs found in the ATAC assay.")
}

# Step 11. Display heatmap
draw(ht)
