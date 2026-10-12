
# ============================================================
# Differential Histone Modification Peak Analysis in Rice
# H3K4me3 and H3K27ac: Heat vs. Control Across Cell Types
# ============================================================

# Step 1. Load required libraries
library(Seurat)
library(Signac)
library(dplyr)
library(tibble)

# Step 2. Set input and output directories
data_dir <- "data"
output_dir <- "marker_results"

dir.create(
  file.path(output_dir, "H3K4me3"),
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  file.path(output_dir, "H3K27ac"),
  showWarnings = FALSE,
  recursive = TRUE
)

# Step 3. Read Seurat objects
h3k4 <- readRDS(
  file.path(data_dir, "h3k4.rds")
)

h3k27 <- readRDS(
  file.path(data_dir, "h3k27.rds")
)

# Step 4. Normalize ATAC assays and define identities
prepare_object <- function(obj) {

  DefaultAssay(obj) <- "ATAC"

  obj <- RunTFIDF(obj)

  obj$celltype <- paste(
    obj$batch,
    obj$cell_type,
    sep = "_"
  )

  Idents(obj) <- "celltype"

  return(obj)
}

h3k4 <- prepare_object(h3k4)
h3k27 <- prepare_object(h3k27)

# Step 5. Define cell types and output labels
cell_info <- data.frame(
  id = c("me", "lp", "gc", "vs", "ep", "fb"),
  cell_type = c(
    "Mesophyll",
    "Large parenchyma",
    "Guard cell",
    "Vascular",
    "Epidermis",
    "Fiber"
  ),
  file_label = c(
    "Mesophyll",
    "LargeParenchyma",
    "GuardCell",
    "Vascular",
    "Epidermis",
    "Fiber"
  )
)

# ============================================================
# Part 1. Differential Peak Analysis
# ============================================================

# Step 6. Define differential analysis function
run_da <- function(obj, mark, prefix, output_subdir) {

  da_list <- lapply(seq_len(nrow(cell_info)), function(i) {

    ct <- cell_info$cell_type[i]

    heat_group <- paste0("heat_", mark, "_", ct)
    ctrl_group <- paste0("ctrl_", mark, "_", ct)

    # Check whether both groups exist
    if (!all(c(heat_group, ctrl_group) %in% levels(Idents(obj)))) {
      stop(paste("Missing cell groups:", heat_group, ctrl_group))
    }

    # Differential peak analysis using logistic regression
    da <- FindMarkers(
      object = obj,
      ident.1 = heat_group,
      ident.2 = ctrl_group,
      assay = "ATAC",
      test.use = "LR",
      min.pct = 0.1,
      latent.vars = c("nCount_ATAC", "nFeature_ATAC"),
      logfc.threshold = 0.25,
      only.pos = FALSE,
      slot = "data"
    )

    # Save differential peak results
    write.csv(
      da,
      file = file.path(
        output_dir,
        output_subdir,
        paste0(prefix, "_", cell_info$file_label[i], ".csv")
      ),
      row.names = TRUE
    )

    return(da)
  })

  names(da_list) <- cell_info$id

  return(da_list)
}

# Step 7. Identify H3K4me3 differential peaks
k4_list <- run_da(
  obj = h3k4,
  mark = "h3k4me3",
  prefix = "k4",
  output_subdir = "H3K4me3"
)

# Step 8. Identify H3K27ac differential peaks
k27_list <- run_da(
  obj = h3k27,
  mark = "h3k27ac",
  prefix = "k27",
  output_subdir = "H3K27ac"
)

# ============================================================
# Part 2. Filter Significant H3K27ac Differential Peaks
# ============================================================

# Step 9. Filter significant H3K27ac peaks
k27.sig <- lapply(k27_list, function(df) {

  df %>%
    rownames_to_column("peak") %>%
    filter(
      abs(avg_log2FC) > 1,
      p_val_adj < 0.05,
      pct.1 >= 0.1
    )
})

# Step 10. Combine significant peaks across cell types
k27.sig.df <- bind_rows(
  lapply(seq_len(nrow(cell_info)), function(i) {

    ct <- cell_info$id[i]

    k27.sig[[ct]] %>%
      mutate(
        celltype = cell_info$cell_type[i]
      )
  })
)

# Step 11. Save significant H3K27ac peaks
write.csv(
  k27.sig.df,
  file = file.path(
    output_dir,
    "H3K27ac",
    "H3K27ac_significant_differential_peaks.csv"
  ),
  row.names = FALSE
)

# ============================================================
# Part 3. Summarize Differential Peaks
# ============================================================

# Step 12. Count differential peaks per cell type
peak_summary <- data.frame(
  cell_type = cell_info$cell_type,
  H3K4me3 = sapply(k4_list, nrow),
  H3K27ac = sapply(k27_list, nrow),
  H3K27ac_significant = sapply(k27.sig, nrow)
)

# Step 13. Save summary statistics
write.csv(
  peak_summary,
  file = file.path(output_dir, "differential_peak_summary.csv"),
  row.names = FALSE
)

# Step 14. Display results
print(peak_summary)
