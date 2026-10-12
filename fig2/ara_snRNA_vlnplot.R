
# ============================================================
# Marker Gene Expression Violin Plot in Arabidopsis
# Visualization of Cell-Type-Specific Marker Genes
# ============================================================

# Step 1. Load required libraries
library(Seurat)
library(ggplot2)

# Step 2. Read Seurat object
data_dir <- "data"
output_dir <- "results"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

Araun <- readRDS(file.path(data_dir, "Araun.rds"))

# Step 3. Define cell type order
cell_type_order <- c(
  "mosophyll",
  "phloem",
  "myrosin_idioblast",
  "late_S_phase/G2/M_phase",
  "G2/M_phase",
  "late_S_phase",
  "vascular",
  "epidermis",
  "guard_cell",
  "companion_cell",
  "xylem"
)

Araun$cell_type <- factor(
  Araun$cell_type,
  levels = cell_type_order
)

# Step 4. Define marker genes
marker_genes <- c(
  "AT3G63140", "AT2G05100", "AT4G13770", "AT5G25980",
  "AT5G26000", "AT1G18370", "AT4G23800", "AT4G11080",
  "AT5G54640", "AT4G20270", "AT5G43760", "AT1G25450",
  "AT2G21140", "AT5G57350", "AT1G79430", "AT5G19530",
  "AT4G32880", "AT5G60490"
)

# Reverse marker gene order
marker_genes_rev <- rev(marker_genes)

# Step 5. Define cell type colors
cell_colors <- c(
  "#4b6aa8", "#3ca0cf", "#c376a7", "#ad98c3",
  "#cea5c7", "#53738c", "#a5a9b0", "#a78982",
  "#696a6c", "#92699e", "#d69971", "#df5734",
  "#6c408e", "#ac6894", "#d4c2db", "#537eb7",
  "#83ab8e", "#ece399", "#405993"
)

# Step 6. Generate stacked violin plot
p <- VlnPlot(
  Araun,
  assay = "RNA",
  features = marker_genes_rev,
  group.by = "cell_type",
  flip = TRUE,
  stack = TRUE,
  cols = cell_colors
) +
  NoLegend() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 3)
  )

# Step 7. Display plot
print(p)

# Step 8. Save plot
ggsave(
  filename = file.path(output_dir, "Arabidopsis_marker_genes_violin.pdf"),
  plot = p,
  device = "pdf",
  width = 12,
  height = 8,
  units = "in"
)
