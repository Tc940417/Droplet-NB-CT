
# ============================================================
# Rice Single-Cell RNA-seq Analysis
# Seurat Processing, CCA Integration, Clustering, and UMAP
# ============================================================

# Step 1. Load required libraries
library(Seurat)
library(ggplot2)

# Step 2. Set input and output directories
data_dir <- "data/3_clean_mtx"
output_dir <- "results"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Step 3. Read 10X gene expression matrices
sample_info <- data.frame(
  sample = c("ctrl1", "ctrl2", "case1", "case2"),
  group = c("ctrl", "ctrl", "case", "case"),
  folder = c(
    "NIP_TR197_none",
    "NIP_TR198_none",
    "NIP_TR199_45C90min",
    "NIP_TR200_45C90min"
  )
)

# Step 4. Create Seurat objects
seurat_list <- lapply(seq_len(nrow(sample_info)), function(i) {

  counts <- Read10X(
    data.dir = file.path(data_dir, sample_info$folder[i]),
    gene.column = 1
  )

  obj <- CreateSeuratObject(
    counts = counts,
    min.cells = 1,
    min.features = 1,
    project = sample_info$sample[i]
  )

  obj$spl <- sample_info$sample[i]
  obj$group <- sample_info$group[i]

  return(obj)
})

names(seurat_list) <- sample_info$sample

# Step 5. Merge samples
obj <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list)
)

# Step 6. Filter cells based on detected genes
obj <- subset(
  obj,
  subset = nFeature_RNA > 500 & nFeature_RNA < 8000
)

# Step 7. Normalize and scale data
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)

# Step 8. Perform PCA
obj <- RunPCA(obj, npcs = 30)

# Step 9. Perform CCA integration
obj <- IntegrateLayers(
  object = obj,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca"
)

# Step 10. Cluster cells using integrated dimensions
obj <- FindNeighbors(
  obj,
  reduction = "integrated.cca",
  dims = 1:30
)

obj <- FindClusters(
  obj,
  resolution = 0.25
)

# Step 11. Run UMAP
set.seed(123)

obj <- RunUMAP(
  obj,
  reduction = "integrated.cca",
  dims = 1:30
)

# Step 12. Save Seurat object
saveRDS(
  obj,
  file = file.path(
    output_dir,
    "nip.gene_500-8000.resolution_0.25.rds"
  )
)

# Step 13. Visualize UMAP by experimental group
p1 <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "group"
) +
  ggtitle("CCA Integration UMAP - Group")

ggsave(
  filename = file.path(output_dir, "Rice_CCA_UMAP_group.pdf"),
  plot = p1,
  width = 8,
  height = 6
)

# Step 14. Visualize UMAP by sample
p2 <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "spl"
) +
  ggtitle("CCA Integration UMAP - Sample")

ggsave(
  filename = file.path(output_dir, "Rice_CCA_UMAP_sample.pdf"),
  plot = p2,
  width = 8,
  height = 6
)

# Step 15. Visualize UMAP by cluster
p3 <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE
) +
  ggtitle("CCA Integration UMAP - Clusters")

ggsave(
  filename = file.path(output_dir, "Rice_CCA_UMAP_clusters.pdf"),
  plot = p3,
  width = 8,
  height = 6
)

# Step 16. Join RNA layers for differential expression
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj, assay = "RNA")

# Step 17. Identify cluster marker genes
Idents(obj) <- "seurat_clusters"

deg <- FindAllMarkers(
  obj,
  assay = "RNA",
  only.pos = FALSE
)

# Step 18. Save differential expression results
write.csv(
  deg,
  file = file.path(
    output_dir,
    "Rice_cluster_markers_resolution_0.25.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

# Step 19. Save final Seurat object
saveRDS(
  obj,
  file = file.path(
    output_dir,
    "nip.gene_500-8000.resolution_0.25.rds"
  )
)

# Step 20. Display plots
print(p1)
print(p2)
print(p3)
