
# ============================================================
# Gene Activity and Marker Heatmap Analysis in Rice
# Comparison of H3K4me3 and H3K27ac
# ============================================================

# Step 1. Load required libraries
library(Seurat)
library(Signac)
library(dplyr)
library(SCP)

# Step 2. Set input and output directories
data_dir <- "snCT"
output_dir <- "results/gene_activity"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Step 3. Read annotated Seurat objects
sn_H3K4me3 <- readRDS(
  file.path(data_dir, "sn_H3K4me3.anno.rds")
)

sn_H3K27ac <- readRDS(
  file.path(data_dir, "sn_h3k27ac.anno.rds")
)

# Step 4. Calculate gene activity scores
gene.activities1 <- GeneActivity(sn_H3K4me3)
gene.activities2 <- GeneActivity(sn_H3K27ac)

# Step 5. Create RNA assays from gene activity matrices
sn_H3K4me3[["RNA"]] <- CreateAssayObject(
  counts = gene.activities1
)

sn_H3K27ac[["RNA"]] <- CreateAssayObject(
  counts = gene.activities2
)

# Step 6. Normalize gene activity scores
sn_H3K4me3 <- NormalizeData(
  object = sn_H3K4me3,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = median(sn_H3K4me3$nCount_RNA)
)

sn_H3K27ac <- NormalizeData(
  object = sn_H3K27ac,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = median(sn_H3K27ac$nCount_RNA)
)

# Step 7. Set default assays and cell identities
DefaultAssay(sn_H3K4me3) <- "RNA"
DefaultAssay(sn_H3K27ac) <- "RNA"

# Define cell type order
cell_type_order <- c(
  "Mesophyll",
  "Large parenchyma",
  "Epidermis",
  "Vascular",
  "Fiber",
  "Guard cell"
)

plot_order <- rev(cell_type_order)

sn_H3K4me3$cell_type_plot <- factor(
  sn_H3K4me3$cell_type,
  levels = plot_order
)

sn_H3K27ac$cell_type_plot <- factor(
  sn_H3K27ac$cell_type,
  levels = plot_order
)

# Retain only the six selected cell types
sn_H3K4me3 <- subset(
  sn_H3K4me3,
  subset = !is.na(cell_type_plot)
)

sn_H3K27ac <- subset(
  sn_H3K27ac,
  subset = !is.na(cell_type_plot)
)

Idents(sn_H3K4me3) <- "cell_type_plot"
Idents(sn_H3K27ac) <- "cell_type_plot"

# Step 8. Save Seurat objects with gene activity assays
saveRDS(
  sn_H3K4me3,
  file = file.path(output_dir, "sn_H3K4me3.anno.addRNA.rds")
)

saveRDS(
  sn_H3K27ac,
  file = file.path(output_dir, "sn_H3K27ac.anno.addRNA.rds")
)

# ============================================================
# H3K27ac Gene Activity Analysis
# ============================================================

# Step 9. Identify H3K27ac cell-type markers
sn_H3K27ac_markers <- FindAllMarkers(
  sn_H3K27ac,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.15,
  logfc.threshold = 0.15
)

# Select the cell type with the highest log2FC per gene
sn_H3K27ac_unique_markers <- sn_H3K27ac_markers %>%
  group_by(gene) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  mutate(
    cluster = factor(cluster, levels = plot_order)
  ) %>%
  arrange(cluster)

# Step 10. Generate H3K27ac gene activity heatmap
sn_H3K27ac_ht <- GroupHeatmap(
  srt = sn_H3K27ac,
  assay = "RNA",
  slot = "data",
  lib_normalize = FALSE,
  group.by = "cell_type_plot",
  features = sn_H3K27ac_unique_markers$gene,
  features_label = sn_H3K27ac_unique_markers$gene,
  height = 6.5,
  width = 4,
  group_palette = "simspec",
  heatmap_palcolor = c("#64a776", "white", "#537eb7"),
  column_title = "sn_H3K27ac Average expression heatmap"
)

# Step 11. Save H3K27ac results
write.csv(
  sn_H3K27ac_markers,
  file = file.path(output_dir, "sn_H3K27ac_markers.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  sn_H3K27ac_unique_markers,
  file = file.path(output_dir, "sn_H3K27ac_unique_markers.csv"),
  row.names = FALSE,
  quote = FALSE
)

pdf(
  file.path(output_dir, "sn_H3K27ac_gene_activity_heatmap.pdf"),
  width = 8,
  height = 8
)

print(sn_H3K27ac_ht$plot)
dev.off()

# ============================================================
# H3K4me3 Gene Activity Analysis
# ============================================================

# Step 12. Identify H3K4me3 cell-type markers
sn_H3K4me3_markers <- FindAllMarkers(
  sn_H3K4me3,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# Select the cell type with the highest log2FC per gene
sn_H3K4me3_unique_markers <- sn_H3K4me3_markers %>%
  group_by(gene) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  mutate(
    cluster = factor(cluster, levels = plot_order)
  ) %>%
  arrange(cluster)

# Step 13. Generate H3K4me3 gene activity heatmap
sn_H3K4me3_ht <- GroupHeatmap(
  srt = sn_H3K4me3,
  assay = "RNA",
  slot = "data",
  lib_normalize = FALSE,
  group.by = "cell_type_plot",
  features = sn_H3K4me3_unique_markers$gene,
  features_label = sn_H3K4me3_unique_markers$gene,
  height = 6.5,
  width = 4,
  group_palette = "simspec",
  heatmap_palcolor = c("#64a776", "white", "#537eb7"),
  column_title = "sn_H3K4me3 Average expression heatmap"
)

# Step 14. Save H3K4me3 results
write.csv(
  sn_H3K4me3_markers,
  file = file.path(output_dir, "sn_H3K4me3_markers.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  sn_H3K4me3_unique_markers,
  file = file.path(output_dir, "sn_H3K4me3_unique_markers.csv"),
  row.names = FALSE,
  quote = FALSE
)

pdf(
  file.path(output_dir, "sn_H3K4me3_gene_activity_heatmap.pdf"),
  width = 8,
  height = 8
)

print(sn_H3K4me3_ht$plot)
dev.off()

# Step 15. Display heatmaps
print(sn_H3K27ac_ht$plot)
print(sn_H3K4me3_ht$plot)
