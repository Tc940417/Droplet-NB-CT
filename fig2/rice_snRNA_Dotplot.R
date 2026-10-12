
# ============================================================
# Marker Gene Expression DotPlot in Rice
# Visualization of Cell-Type-Specific Marker Genes
# ============================================================

# Step 1. Load required libraries
library(Seurat)
library(ggplot2)

# Step 2. Read annotated Seurat object
data_dir <- "data"
output_dir <- "results"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

nip <- readRDS(file.path(data_dir, "nip_annotated.rds"))

# Step 3. Define cell type order
cell_type_order <- c(
  "Mesophyll",
  "Large parenchyma",
  "Guard cell",
  "Fiber",
  "Epidermis",
  "Vascular"
)

# Step 4. Define marker genes
marker_genes <- c(
  "LOC-Os05g02420",
  "LOC-Os06g47640",
  "LOC-Os04g58800",
  "LOC-Os03g55240",
  "LOC-Os02g41860",
  "LOC-Os04g20164",
  "LOC-Os12g02320",
  "LOC-Os11g10310",
  "LOC-Os11g10320",
  "LOC-Os09g10054",
  "LOC-Os08g43670",
  "LOC-Os07g46852",
  "LOC-Os01g15010",
  "LOC-Os10g42190",
  "LOC-Os01g64730",
  "LOC-Os11g13890",
  "LOC-Os09g17740",
  "LOC-Os08g33820",
  "LOC-Os07g38960"
)

# Step 5. Define color palette
dot_colors <- c(
  "#330066",
  "#336699",
  "#66CC66",
  "#FFCC33"
)

# Step 6. Generate marker gene DotPlot
p <- DotPlot(
  nip,
  assay = "RNA",
  features = marker_genes,
  group.by = "cell_type"
) +
  scale_y_discrete(limits = cell_type_order) +
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

# Step 7. Display plot
print(p)

# Step 8. Save plot
ggsave(
  filename = file.path(output_dir, "Rice_marker_genes_dotplot.pdf"),
  plot = p,
  device = "pdf",
  width = 10,
  height = 6,
  units = "in"
)
