
# ============================================================
# Arabidopsis Single-Cell RNA-seq Analysis
# Seurat Object Creation, Quality Control, and Clustering
# ============================================================

# Step 1. Load required libraries
library(Seurat)

# Step 2. Set input and output directories
data_dir <- "data/3_clean_mtx"
output_dir <- "results"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Step 3. Read 10X gene expression matrices
cnt91 <- Read10X(
  data.dir = file.path(data_dir, "Ara_TR091_none"),
  gene.column = 1
)

cnt92 <- Read10X(
  data.dir = file.path(data_dir, "Ara_TR092_none"),
  gene.column = 1
)

# Step 4. Create Seurat objects
obj91 <- CreateSeuratObject(
  counts = cnt91,
  min.features = 1,
  min.cells = 1,
  project = "ath_91"
)

obj92 <- CreateSeuratObject(
  counts = cnt92,
  min.features = 1,
  min.cells = 1,
  project = "ath_92"
)

# Add sample information
obj91$spl <- "ath_91"
obj92$spl <- "ath_92"

# Step 5. Merge samples
obj <- merge(
  obj91,
  y = obj92,
  add.cell.ids = c("ath_91", "ath_92")
)

# Join RNA assay layers (Seurat v5)
obj <- JoinLayers(obj)

# Step 6. Filter cells based on detected genes
obj <- subset(
  obj,
  subset = nFeature_RNA > 500 & nFeature_RNA < 5000
)

# Step 7. Normalize and scale data
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)

# Step 8. Perform dimensionality reduction
obj <- RunPCA(obj, npcs = 30)

# Step 9. Cluster cells
obj <- FindNeighbors(obj, dims = 1:20)
obj <- FindClusters(obj, resolution = 0.5)

# Step 10. Run UMAP
set.seed(123)
obj <- RunUMAP(obj, dims = 1:20)

# Step 11. Identify cluster marker genes
deg <- FindAllMarkers(
  obj,
  only.pos = FALSE
)

# Step 12. Save results
saveRDS(
  obj,
  file = file.path(output_dir, "Ara.rds")
)

write.csv(
  deg,
  file = file.path(output_dir, "Arabidopsis_cluster_markers.csv"),
  row.names = FALSE
)

# Step 13. Visualize UMAP
p <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE
)

print(p)

ggsave(
  filename = file.path(output_dir, "Arabidopsis_UMAP.pdf"),
  plot = p,
  width = 8,
  height = 6
)
