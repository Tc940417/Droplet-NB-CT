
# ============================================================
# H3K4me3 Differential Peak Heatmap Analysis in Rice
# Cell-Type-Specific Peak Accessibility (Z-score)
# ============================================================

# Step 1. Load required libraries
library(Seurat)
library(Signac)
library(ComplexHeatmap)
library(circlize)
library(grid)

# Step 2. Set input and output directories
data_dir <- "data"
output_dir <- "results/peak_heatmap"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Step 3. Read H3K4me3 Seurat object
h3k4 <- readRDS(
  file.path(data_dir, "ct.h3k4.rds")
)

DefaultAssay(h3k4) <- "ATAC"

# Step 4. Normalize ATAC data using TF-IDF
h3k4 <- RunTFIDF(h3k4)

# Step 5. Read differential peak results
h3k4.da <- read.csv(
  file.path(data_dir, "ctH3K4_FindAllMarkers_LR.csv"),
  header = TRUE,
  stringsAsFactors = FALSE
)

# Step 6. Retain differential peaks present in the ATAC assay
peak_features <- unique(h3k4.da$gene)
peak_features <- intersect(
  peak_features,
  rownames(h3k4[["ATAC"]])
)

if (length(peak_features) == 0) {
  stop("No differential peaks match the ATAC assay features.")
}

# Step 7. Annotate peaks with closest genes
cf <- ClosestFeature(
  h3k4,
  regions = peak_features
)

write.csv(
  cf,
  file = file.path(output_dir, "H3K4me3_peak_annotations.csv"),
  row.names = FALSE
)

# Step 8. Define cell type order
cell_type_order <- c(
  "Epidermis",
  "Guard cell",
  "Mesophyll",
  "Large parenchyma",
  "Vascular",
  "Fiber"
)

# Step 9. Calculate mean peak signals for each cell type
# Use TF-IDF normalized ATAC data directly
atac_matrix <- GetAssayData(
  h3k4,
  assay = "ATAC",
  layer = "data"
)

# Retain selected cell types
cell_types <- as.character(h3k4$cell_type)
names(cell_types) <- colnames(h3k4)

valid_cells <- names(cell_types)[
  !is.na(cell_types) & cell_types %in% cell_type_order
]

atac_matrix <- atac_matrix[
  peak_features,
  valid_cells,
  drop = FALSE
]

cell_types <- cell_types[valid_cells]

# Calculate mean signal per peak per cell type
peak_avg_matrix <- sapply(
  cell_type_order,
  function(ct) {
    cells <- names(cell_types)[cell_types == ct]
    if (length(cells) == 0) {
      stop(paste("No cells found for:", ct))
    }
    Matrix::rowMeans(atac_matrix[, cells, drop = FALSE])
  }
)

peak_avg_matrix <- as.matrix(peak_avg_matrix)

# Step 10. Calculate Z-scores across cell types
# Rows: cell types; columns: peaks
peak_zscore <- t(scale(t(peak_avg_matrix)))

# Remove peaks with zero variance or invalid Z-scores
peak_zscore <- peak_zscore[
  rowSums(!is.finite(peak_zscore)) == 0,
  ,
  drop = FALSE
]

if (nrow(peak_zscore) == 0) {
  stop("No peaks remain after Z-score filtering.")
}

sig.region.mean.scale <- t(peak_zscore)

# Step 11. Set cell type order
sig.region.mean.scale <- sig.region.mean.scale[
  cell_type_order,
  ,
  drop = FALSE
]

# Step 12. Order peaks by the cell type with maximum signal
dominant_group <- apply(
  sig.region.mean.scale,
  2,
  which.max
)

col_order <- order(dominant_group)

sig.region.mean.scale <- sig.region.mean.scale[
  ,
  col_order,
  drop = FALSE
]

# Step 13. Define peak groups
col_group <- apply(
  sig.region.mean.scale,
  2,
  function(x) {
    rownames(sig.region.mean.scale)[which.max(x)]
  }
)

col_group <- factor(
  col_group,
  levels = cell_type_order
)

# Step 14. Define heatmap color palette
lim <- 2

col_fun <- colorRamp2(
  c(-lim, -1, 0, 1, lim),
  c(
    "#6F8F7A",
    "#B8CFC1",
    "#FFFFFF",
    "#B7C9D9",
    "#5B7FA6"
  )
)

# Step 15. Generate H3K4me3 peak heatmap
ht <- Heatmap(
  sig.region.mean.scale,
  name = "Z-score",
  col = col_fun,

  cluster_rows = FALSE,
  cluster_columns = FALSE,

  column_split = col_group,
  column_title = NULL,

  show_column_names = FALSE,

  row_names_gp = gpar(
    fontsize = 8,
    fontface = "bold"
  ),

  width = unit(12, "cm"),
  height = unit(8, "cm"),

  border = TRUE,

  heatmap_legend_param = list(
    title = "Z-score",
    at = c(-2, -1, 0, 1, 2),
    labels = c("-2", "-1", "0", "1", "2"),
    title_gp = gpar(fontsize = 10),
    labels_gp = gpar(fontsize = 9),
    legend_height = unit(4, "cm")
  )
)

# Step 16. Save heatmap
pdf(
  file.path(output_dir, "Rice_H3K4me3_DA_peak_heatmap.pdf"),
  width = 10,
  height = 6
)

draw(ht, heatmap_legend_side = "right")
dev.off()

# Step 17. Save Z-score matrix
write.csv(
  sig.region.mean.scale,
  file = file.path(output_dir, "Rice_H3K4me3_peak_zscore.csv"),
  row.names = TRUE,
  quote = FALSE
)

# Step 18. Display heatmap
draw(ht, heatmap_legend_side = "right")
