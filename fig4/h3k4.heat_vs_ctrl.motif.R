
# ============================================================
# Differential Motif Activity Analysis in Rice
# H3K4me3 chromVAR: Heat vs. Control Across Cell Types
# ============================================================

# Step 1. Load required libraries
library(Seurat)
library(Signac)
library(dplyr)
library(tibble)

# Step 2. Set input and output directories
data_dir <- "data"
output_dir <- "results/H3K4me3_differential_motifs"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Step 3. Read H3K4me3 Seurat object
h3k4 <- readRDS(
  file.path(data_dir, "h3k4.rds")
)

# Step 4. Set chromVAR assay and cell identities
DefaultAssay(h3k4) <- "chromvar"

h3k4$celltype <- paste(
  h3k4$batch,
  h3k4$cell_type,
  sep = "-"
)

Idents(h3k4) <- "celltype"

# Step 5. Define cell types
cell_types <- c(
  "Mesophyll",
  "Large parenchyma",
  "Epidermis",
  "Vascular",
  "Fiber",
  "Guard cell"
)

# ============================================================
# Differential Motif Activity Analysis
# ============================================================

# Step 6. Compare heat and control samples for each cell type
da.motif_list <- lapply(cell_types, function(ct) {

  heat_group <- paste0("heat_h3k4me3-", ct)
  ctrl_group <- paste0("ctrl_h3k4me3-", ct)

  FindMarkers(
    object = h3k4,
    ident.1 = heat_group,
    ident.2 = ctrl_group,
    assay = "chromvar",
    slot = "data",
    only.pos = FALSE,
    mean.fxn = rowMeans,
    fc.name = "avg_diff"
  )
})

names(da.motif_list) <- cell_types

# Step 7. Combine differential motif results
da.motif <- bind_rows(
  lapply(names(da.motif_list), function(ct) {

    df <- da.motif_list[[ct]]

    if (nrow(df) == 0) {
      return(NULL)
    }

    df <- df %>%
      rownames_to_column("motif") %>%
      mutate(
        cell_type = ct,
        padj_plot = pmax(p_val_adj, 1e-303),
        neglog10_padj = -log10(padj_plot)
      )

    # Calculate detection fraction difference if available
    if (all(c("pct.1", "pct.2") %in% colnames(df))) {
      df$pct_diff <- abs(df$pct.1 - df$pct.2)
    }

    return(df)
  })
)

# Step 8. Save differential motif results for each cell type
for (ct in names(da.motif_list)) {

  file_name <- paste0(
    "H3K4me3_heat_vs_ctrl_",
    gsub(" ", "_", ct),
    "_differential_motif_activity.csv"
  )

  write.csv(
    da.motif_list[[ct]],
    file = file.path(output_dir, file_name),
    row.names = TRUE,
    quote = FALSE
  )
}

# Step 9. Save combined differential motif results
write.csv(
  da.motif,
  file = file.path(
    output_dir,
    "H3K4me3_all_celltypes_differential_motifs.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

# Step 10. Display results
print(head(da.motif))
