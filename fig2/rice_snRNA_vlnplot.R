
# ============================================================
# Marker Gene Expression Violin Plot in Rice
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

# Step 4. Remove unknown cells and set cell type order
nip_filtered <- subset(
  nip,
  subset = cell_type %in% cell_type_order
)

nip_filtered$cell_type <- factor(
  nip_filtered$cell_type,
  levels = cell_type_order
)

# Step 5. Define marker genes
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

# Step 6. Define cell type colors
cell_colors <- c(
  "Mesophyll" = "#927c9a",
  "Large parenchyma" = "#3674a2",
  "Guard cell" = "#9f8d89",
  "Fiber" = "#72567a",
  "Epidermis" = "#63a3b8",
  "Vascular" = "#c4daec"
)

# Step 7. Generate stacked violin plot
p <- VlnPlot(
  nip_filtered,
  assay = "RNA",
  features = marker_genes,
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

# Step 8. Display plot
print(p)

# Step 9. Save plot
ggsave(
  filename = file.path(output_dir, "Rice_marker_genes_violin.pdf"),
  plot = p,
  device = "pdf",
  width = 12,
  height = 8,
  units = "in"
)
