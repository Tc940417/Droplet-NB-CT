
# ============================================================
# Differential H3K4me3 Peak Analysis in Rice
# Single-Cell Logistic Regression
# Heat vs. Control in Vascular Cells
# ============================================================

# Step 1. Load required libraries
library(Seurat)
library(Signac)

# Step 2. Set input and output directories
data_dir <- "data"
output_dir <- "results/H3K4me3_differential_peaks"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Step 3. Read H3K4me3 Seurat object
h3k4 <- readRDS(
  file.path(data_dir, "h3k4.rds")
)

DefaultAssay(h3k4) <- "ATAC"

# Step 4. Normalize ATAC data using TF-IDF
h3k4 <- RunTFIDF(h3k4)

# Step 5. Define experimental conditions
h3k4$condition <- as.character(h3k4$batch)

h3k4$condition[
  h3k4$condition == "ctrl_h3k4me3"
] <- "ctrl"

h3k4$condition[
  h3k4$condition == "heat_h3k4me3"
] <- "heat"

# Define cell type and condition identities
h3k4$label <- paste(
  h3k4$cell_type,
  h3k4$condition,
  sep = "-"
)

Idents(h3k4) <- "label"

# Step 6. Select target cell type
ct <- "Vascular"

# ============================================================
# Part 1. Pseudobulk Differential Peak Analysis (DESeq2)
# ============================================================

# Step 7. Subset target cell type
h3k4_ct <- subset(
  h3k4,
  subset = cell_type == ct
)

# Step 8. Prepare biological replicate information
# dataset must identify independent biological samples
meta <- h3k4_ct@meta.data

sample_info <- unique(
  meta[, c("dataset", "condition"), drop = FALSE]
)

if (anyNA(sample_info)) {
  stop("Missing dataset or condition information.")
}

if (!all(sample_info$condition %in% c("ctrl", "heat"))) {
  stop("Unexpected condition labels.")
}

# Assign a unique pseudobulk ID to each sample
sample_info$pb_id <- paste0(
  "sample_", seq_len(nrow(sample_info))
)

sample_key <- paste(
  sample_info$dataset,
  sample_info$condition,
  sep = "__"
)

cell_key <- paste(
  meta$dataset,
  meta$condition,
  sep = "__"
)

h3k4_ct$pb_id <- sample_info$pb_id[
  match(cell_key, sample_key)
]

# Step 9. Aggregate raw peak counts by biological replicate
pb_counts <- AggregateExpression(
  object = h3k4_ct,
  assays = "ATAC",
  group.by = "pb_id",
  return.seurat = FALSE
)$ATAC

# Step 10. Prepare pseudobulk metadata
coldata <- sample_info[
  match(colnames(pb_counts), sample_info$pb_id),
  ,
  drop = FALSE
]

rownames(coldata) <- colnames(pb_counts)

coldata$condition <- factor(
  coldata$condition,
  levels = c("ctrl", "heat")
)

# Check biological replicates
sample_counts <- table(coldata$condition)
print(sample_counts)

if (any(sample_counts < 2)) {
  stop("At least two biological replicates per condition are required.")
}

# Step 11. Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(pb_counts)),
  colData = coldata,
  design = ~ condition
)

# Step 12. Filter low-count peaks
dds <- dds[
  rowSums(counts(dds)) > 10,
]

# Step 13. Perform differential peak analysis
dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c("condition", "heat", "ctrl")
)

# Step 14. Save DESeq2 results
res_df <- as.data.frame(res)
res_df$peak <- rownames(res_df)

res_df <- res_df[
  , c("peak", setdiff(colnames(res_df), "peak"))
]

write.csv(
  res_df,
  file = file.path(
    output_dir,
    paste0("H3K4me3_", ct, "_DESeq2.csv")
  ),
  row.names = FALSE
)

write.csv(
  coldata,
  file = file.path(
    output_dir,
    paste0("H3K4me3_", ct, "_pseudobulk_samples.csv")
  ),
  row.names = FALSE
)

# ============================================================
# Part 2. Single-Cell Differential Peak Analysis (LR)
# ============================================================

# Step 15. Identify differential peaks using logistic regression
da.lr <- FindMarkers(
  object = h3k4,
  ident.1 = paste0(ct, "-heat"),
  ident.2 = paste0(ct, "-ctrl"),
  assay = "ATAC",
  test.use = "LR",
  latent.vars = c("nCount_ATAC", "nFeature_ATAC"),
  min.pct = 0.05,
  logfc.threshold = 0,
  only.pos = FALSE,
  slot = "data"
)

# Step 16. Save LR differential peak results
da.lr_df <- data.frame(
  peak = rownames(da.lr),
  da.lr,
  row.names = NULL,
  check.names = FALSE
)

write.csv(
  da.lr_df,
  file = file.path(
    output_dir,
    paste0("H3K4me3_", ct, "_LR.csv")
  ),
  row.names = FALSE
)

# Step 17. Display results
print(head(res_df))
print(head(da.lr_df))
