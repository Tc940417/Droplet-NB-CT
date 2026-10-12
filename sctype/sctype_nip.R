
# ============================================================
# Automated Cell Type Annotation in Rice
# scType Annotation Using a Custom Marker Gene Database
# ============================================================

# Step 1. Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(HGNChelper)
library(openxlsx)

# Step 2. Load scType functions
source("scripts/gene_sets_prepare.R")
source("scripts/sctype_score.R")

# Step 3. Set input and output directories
data_dir <- "data"
output_dir <- "results/scType"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Step 4. Read Seurat object
nip.sctype <- readRDS(
  file.path(data_dir, "nip.gene_500-5000.resolution_0.25.rds")
)

# Step 5. Extract scaled RNA expression matrix
DefaultAssay(nip.sctype) <- "RNA"

scaled_matrix <- LayerData(
  nip.sctype,
  assay = "RNA",
  layer = "scale.data"
)

# Verify that scaled data include all cells
if (!all(colnames(nip.sctype) %in% colnames(scaled_matrix))) {
  stop("The RNA scale.data layer does not contain all cells.")
}

scaled_matrix <- scaled_matrix[
  ,
  colnames(nip.sctype),
  drop = FALSE
]

# Step 6. Load custom marker gene database
marker_db <- file.path(
  data_dir,
  "ricemerge_DB -.xlsx"
)

tissue <- "Leaf"

gs_list <- gene_sets_prepare(
  marker_db,
  tissue
)

# Step 7. Calculate scType scores
es.max <- sctype_score(
  scRNAseqData = as.matrix(scaled_matrix),
  scaled = TRUE,
  gs = gs_list$gs_positive,
  gs2 = NULL,
  gene_names_to_uppercase = FALSE
)

# ============================================================
# Part 1. Calculate Cluster-Level scType Scores
# ============================================================

# Step 8. Extract Seurat cluster information
clusters <- as.character(nip.sctype$seurat_clusters)
names(clusters) <- colnames(nip.sctype)

# Step 9. Calculate cell type scores for each cluster
cluster_results <- lapply(unique(clusters), function(cl) {

  cells <- names(clusters)[clusters == cl]

  # Sum scType scores across cells in each cluster
  cluster_scores <- sort(
    rowSums(es.max[, cells, drop = FALSE]),
    decreasing = TRUE
  )

  # Retain top 12 candidate cell types
  head(
    data.frame(
      cluster = cl,
      type = names(cluster_scores),
      scores = as.numeric(cluster_scores),
      ncells = length(cells)
    ),
    12
  )
})

cluster_results <- bind_rows(cluster_results)

# ============================================================
# Part 2. Assign Cell Types to Clusters
# ============================================================

# Step 10. Select the highest-scoring cell type per cluster
sctype_scores <- cluster_results %>%
  group_by(cluster) %>%
  slice_max(
    order_by = scores,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()

# Step 11. Assign low-confidence clusters as Unknown
sctype_scores <- sctype_scores %>%
  mutate(
    type = ifelse(
      scores < ncells / 4,
      "Unknown",
      type
    )
  )

# Step 12. Display annotation results
print(
  sctype_scores[, c("cluster", "type", "scores", "ncells")]
)

# ============================================================
# Part 3. Add Cell Type Annotations to Seurat Object
# ============================================================

# Step 13. Map cluster annotations to individual cells
cluster_annotation <- setNames(
  sctype_scores$type,
  sctype_scores$cluster
)

nip.sctype$customclassif <- unname(
  cluster_annotation[as.character(nip.sctype$seurat_clusters)]
)

# ============================================================
# Part 4. Visualize Cell Type Annotations
# ============================================================

# Step 14. Generate annotated UMAP
p <- DimPlot(
  nip.sctype,
  reduction = "umap",
  group.by = "customclassif",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("scType Cell Type Annotation in Rice")

# Step 15. Save UMAP
ggsave(
  filename = file.path(output_dir, "Rice_scType_UMAP.pdf"),
  plot = p,
  device = "pdf",
  width = 8,
  height = 6
)

# ============================================================
# Part 5. Save Annotation Results
# ============================================================

# Step 16. Save annotated Seurat object
saveRDS(
  nip.sctype,
  file = file.path(output_dir, "nip.sctype.anno.rds")
)

# Step 17. Save cluster annotation results
write.csv(
  sctype_scores,
  file = file.path(output_dir, "Rice_scType_cluster_annotations.csv"),
  row.names = FALSE
)

# Step 18. Save all cluster-level scores
write.csv(
  cluster_results,
  file = file.path(output_dir, "Rice_scType_cluster_scores.csv"),
  row.names = FALSE
)

# Step 19. Display UMAP
print(p)
