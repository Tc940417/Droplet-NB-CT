
# ============================================================
# Differential Gene Expression Analysis in Rice
# MAST: Heat vs. Control Across Six Cell Types
# ============================================================

# Step 1. Load required libraries
library(Seurat)
library(MAST)

# Step 2. Set input and output directories
data_dir <- "data"
output_dir <- "results/rice_snRNA_DEG"
split_dir <- file.path(output_dir, "DEGs_split")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(split_dir, showWarnings = FALSE, recursive = TRUE)

# Step 3. Read annotated scRNA-seq Seurat object
nip <- readRDS(
  file.path(data_dir, "nip_annotated.rds")
)

# Step 4. Prepare RNA assay for differential analysis
DefaultAssay(nip) <- "RNA"
nip <- JoinLayers(nip, assay = "RNA")

# Define cell type and treatment identities
nip$label <- paste(
  nip$cell_type,
  nip$group,
  sep = "-"
)

Idents(nip) <- "label"

# Step 5. Define cell types
cell_types <- c(
  Mesophyll = "Mesophyll",
  LargeParenchyma = "Large parenchyma",
  Epidermis = "Epidermis",
  Vascular = "Vascular",
  GuardCell = "Guard cell",
  Fiber = "Fiber"
)

# ============================================================
# Part 1. Differential Gene Expression Analysis
# ============================================================

# Step 6. Perform MAST analysis for each cell type
deg_list <- lapply(names(cell_types), function(ct) {

  cell_type <- cell_types[[ct]]

  case_group <- paste0(cell_type, "-case")
  ctrl_group <- paste0(cell_type, "-ctrl")

  # Check whether both groups exist
  if (!all(c(case_group, ctrl_group) %in% levels(Idents(nip)))) {
    stop(paste("Missing cell groups:", case_group, ctrl_group))
  }

  # Identify differentially expressed genes
  deg <- FindMarkers(
    object = nip,
    ident.1 = case_group,
    ident.2 = ctrl_group,
    assay = "RNA",
    test.use = "MAST",
    min.pct = 0.1,
    logfc.threshold = 0.25,
    latent.vars = "nCount_RNA",
    slot = "data",
    only.pos = FALSE
  )

  # Save complete DEG results
  write.csv(
    deg,
    file = file.path(
      output_dir,
      paste0(ct, "_case_vs_ctrl.csv")
    ),
    row.names = TRUE
  )

  return(deg)
})

names(deg_list) <- names(cell_types)

# ============================================================
# Part 2. Classify Upregulated and Downregulated Genes
# ============================================================

# Step 7. Define DEG filtering thresholds
padj_cutoff <- 0.05
log2FC_cutoff <- 1

# Step 8. Classify DEGs
deg_split <- lapply(deg_list, function(df) {

  significant <- is.finite(df$p_val_adj) &
    is.finite(df$avg_log2FC) &
    df$p_val_adj < padj_cutoff

  up <- significant & df$avg_log2FC >= log2FC_cutoff
  down <- significant & df$avg_log2FC <= -log2FC_cutoff

  list(
    up = df[up, , drop = FALSE],
    down = df[down, , drop = FALSE],
    ns = df[!(up | down), , drop = FALSE]
  )
})

# Step 9. Save classified DEGs for each cell type
for (ct in names(deg_split)) {

  ct_dir <- file.path(split_dir, ct)

  dir.create(
    ct_dir,
    showWarnings = FALSE,
    recursive = TRUE
  )

  for (category in c("up", "down", "ns")) {

    write.csv(
      deg_split[[ct]][[category]],
      file = file.path(ct_dir, paste0(category, ".csv")),
      row.names = TRUE
    )
  }
}

# ============================================================
# Part 3. Summarize Differential Gene Expression
# ============================================================

# Step 10. Count DEGs for each cell type
deg_summary <- data.frame(
  cell_type = names(deg_split),
  n_tested = sapply(deg_list, nrow),
  up = sapply(deg_split, function(x) nrow(x$up)),
  down = sapply(deg_split, function(x) nrow(x$down)),
  ns = sapply(deg_split, function(x) nrow(x$ns)),
  row.names = NULL
)

# Step 11. Save DEG summary
write.csv(
  deg_summary,
  file = file.path(output_dir, "Rice_DEG_summary.csv"),
  row.names = FALSE
)

# Step 12. Save complete DEG results as an RDS file
saveRDS(
  deg_list,
  file = file.path(output_dir, "Rice_DEG_results_MAST.rds")
)

# Step 13. Display summary
print(deg_summary)
