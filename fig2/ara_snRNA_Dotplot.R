
# ============================================================
# Marker Gene Expression DotPlot in Arabidopsis
# Visualization of Cell-Type-Specific Marker Genes
# ============================================================
Ara <- readRDS("data/Ara.rds")
# Step 1. Load required libraries
library(Seurat)
library(ggplot2)

# Step 2. Define cell type order
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

# Step 3. Define marker genes
marker_genes <- c(
  "AT3G63140", "AT2G05100", "AT4G13770", "AT5G25980",
  "AT5G26000", "AT1G18370", "AT4G23800", "AT4G11080",
  "AT5G54640", "AT4G20270", "AT5G43760", "AT1G25450",
  "AT2G21140", "AT5G57350", "AT1G79430", "AT5G19530",
  "AT4G32880", "AT5G60490"
)

# Step 4. Define color palette
dot_colors <- c(
  "#330066",
  "#336699",
  "#66CC66",
  "#FFCC33"
)

# Step 5. Generate marker gene DotPlot
p <- DotPlot(
  Ara,
  assay = "RNA",
  features = marker_genes,
  group.by = "cell_type"
) +
  scale_y_discrete(limits = rev(cell_type_order)) +
  coord_flip() +
  scale_color_gradientn(colours = dot_colors) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme(
    axis.text.x = element_text(
      size = 8,
      angle = 45,
      hjust = 1,
      face = "italic"
    ),
    axis.text.y = element_text(size = 8)
  )

# Step 6. Display plot
print(p)

# Step 7. Save plot
dir.create("results", showWarnings = FALSE, recursive = TRUE)

ggsave(
  filename = "results/Arabidopsis_marker_genes_dotplot.pdf",
  plot = p,
  device = "pdf",
  width = 10,
  height = 6,
  units = "in"
)
