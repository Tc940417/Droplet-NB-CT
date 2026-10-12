
# ============================================================
# H3K27ac and RNA Expression Correlation Analysis in Rice
# Cell-Type-Specific Differential Peak and DEG Integration
# ============================================================

# Step 1. Load required libraries
library(Seurat)
library(Signac)
library(dplyr)
library(tibble)
library(ggplot2)

# Step 2. Set input and output directories
rna_dir <- "data/rice_snRNA_DEG"
peak_dir <- "data/H3K27ac_DA"
output_dir <- "results/H3K27ac_RNA_correlation"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Step 3. Define cell types and input files
celltypes <- c("me", "lp", "gc", "vs", "ep", "fb")

rna_files <- c(
  me = "Mesophyll_case_vs_ctrl.csv",
  lp = "LargeParenchyma_case_vs_ctrl.csv",
  gc = "GuardCell_case_vs_ctrl.csv",
  vs = "Vascular_case_vs_ctrl.csv",
  ep = "Epidermis_case_vs_ctrl.csv",
  fb = "Fiber_case_vs_ctrl.csv"
)

# Differential peak results from H3K27ac analysis
peak_files <- c(
  me = "k27.me.csv",
  lp = "k27.lp.csv",
  gc = "k27.gc.csv",
  vs = "k27.vs.csv",
  ep = "k27.ep.csv",
  fb = "k27.fb.csv"
)

# Step 4. Read differential expression results
rna_list <- lapply(rna_files, function(file) {
  read.csv(
    file.path(rna_dir, file),
    row.names = 1,
    check.names = FALSE
  )
})

# Step 5. Read differential peak results
peak_list <- lapply(peak_files, function(file) {
  read.csv(
    file.path(peak_dir, file),
    row.names = 1,
    check.names = FALSE
  )
})

# Step 6. Read H3K27ac Seurat object
h3k27 <- readRDS(
  "data/ct.h3k27.rds"
)

DefaultAssay(h3k27) <- "ATAC"

# ============================================================
# Part 1. Filter Differentially Expressed Genes
# ============================================================

# Step 7. Filter significant RNA DEGs
rna_f_list <- lapply(rna_list, function(df) {
  df %>%
    filter(
      p_val_adj < 0.01,
      abs(avg_log2FC) > 1.5
    ) %>%
    rownames_to_column("gene")
})

# Step 8. Summarize DEG counts
rna_step1_summary <- data.frame(
  celltype = names(rna_f_list),
  n_deg = sapply(rna_f_list, nrow)
)

# ============================================================
# Part 2. Filter Differential H3K27ac Peaks
# ============================================================

# Step 9. Filter significant differential peaks
peak_f_list <- lapply(peak_list, function(df) {
  df %>%
    filter(
      p_val < 0.005,
      abs(avg_log2FC) > 1
    ) %>%
    rownames_to_column("peak")
})

# Step 10. Summarize differential peak counts
peak_step2_summary <- data.frame(
  celltype = names(peak_f_list),
  n_diffpeak = sapply(peak_f_list, nrow)
)

# ============================================================
# Part 3. Annotate Differential Peaks to Nearby Genes
# ============================================================

# Step 11. Annotate peaks with closest genes
peak_anno_list <- lapply(peak_f_list, function(df) {

  if (nrow(df) == 0) {
    df$gene <- character(0)
    df$distance <- numeric(0)
    return(df)
  }

  peak.gr <- StringToGRanges(df$peak)

  cf <- ClosestFeature(
    object = h3k27,
    regions = peak.gr
  )

  # Match annotations using peak coordinates
  cf <- cf %>%
    select(query_region, gene_name, distance) %>%
    distinct(query_region, .keep_all = TRUE)

  df <- df %>%
    left_join(
      cf,
      by = c("peak" = "query_region")
    ) %>%
    rename(
      gene = gene_name
    ) %>%
    mutate(
      gene = gsub("_", "-", gene)
    )

  return(df)
})

# Step 12. Summarize annotated peaks
peak_step3_summary <- data.frame(
  celltype = names(peak_anno_list),
  n_peak_total = sapply(peak_anno_list, nrow),
  n_peak_with_gene = sapply(
    peak_anno_list,
    function(df) sum(!is.na(df$gene))
  )
)

# ============================================================
# Part 4. Calculate Gene-Level H3K27ac Changes
# ============================================================

# Step 13. Aggregate differential peaks by associated genes
peak_gene_list <- lapply(peak_anno_list, function(df) {
  df %>%
    filter(!is.na(gene)) %>%
    group_by(gene) %>%
    summarise(
      peak_logFC = mean(avg_log2FC, na.rm = TRUE),
      n_peak = n(),
      mean_distance = mean(distance, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(abs(peak_logFC) > 1)
})

# Step 14. Summarize peak-associated genes
peak_step4_summary <- data.frame(
  celltype = names(peak_gene_list),
  n_gene_from_peaks = sapply(peak_gene_list, nrow)
)

# ============================================================
# Part 5. Integrate RNA DEGs and H3K27ac Peak Changes
# ============================================================

# Step 15. Match genes between RNA and H3K27ac datasets
merged_list <- lapply(celltypes, function(ct) {

  rna_df <- rna_f_list[[ct]]
  peak_df <- peak_gene_list[[ct]]

  inner_join(
    rna_df %>%
      select(gene, rna_logFC = avg_log2FC),
    peak_df %>%
      select(gene, peak_logFC),
    by = "gene"
  )
})

names(merged_list) <- celltypes

# Step 16. Summarize overlapping genes
merge_step5_summary <- data.frame(
  celltype = names(merged_list),
  n_overlap_gene = sapply(merged_list, nrow)
)

# ============================================================
# Part 6. Spearman Correlation Analysis
# ============================================================

# Step 17. Calculate Spearman correlations
cor_results <- do.call(
  rbind,
  lapply(names(merged_list), function(ct) {

    df <- merged_list[[ct]] %>%
      filter(is.finite(rna_logFC), is.finite(peak_logFC))

    if (
      nrow(df) < 10 ||
      dplyr::n_distinct(df$rna_logFC) < 2 ||
      dplyr::n_distinct(df$peak_logFC) < 2
    ) {
      return(data.frame(
        celltype = ct,
        rho = NA_real_,
        pvalue = NA_real_,
        n_gene = nrow(df)
      ))
    }

    res <- suppressWarnings(cor.test(
      df$rna_logFC,
      df$peak_logFC,
      method = "spearman",
      exact = FALSE
    ))

    data.frame(
      celltype = ct,
      rho = as.numeric(res$estimate),
      pvalue = res$p.value,
      n_gene = nrow(df)
    )
  })
)

cor_results <- cor_results %>%
  arrange(desc(abs(rho)))

# Step 18. Save summary statistics
summary_list <- list(
  RNA_DEGs = rna_step1_summary,
  Differential_Peaks = peak_step2_summary,
  Annotated_Peaks = peak_step3_summary,
  Peak_Associated_Genes = peak_step4_summary,
  Overlapping_Genes = merge_step5_summary,
  Spearman_Correlation = cor_results
)

for (name in names(summary_list)) {
  write.csv(
    summary_list[[name]],
    file = file.path(output_dir, paste0(name, ".csv")),
    row.names = FALSE
  )
}

# Step 19. Save integrated gene-level results
saveRDS(
  merged_list,
  file = file.path(output_dir, "H3K27ac_RNA_merged_genes.rds")
)

for (ct in celltypes) {
  write.csv(
    merged_list[[ct]],
    file = file.path(output_dir, paste0(ct, "_merged_genes.csv")),
    row.names = FALSE
  )
}

# ============================================================
# Part 7. Generate Cell-Type-Specific Correlation Plots
# ============================================================

# Step 20. Define cell type colors
cell_colors <- c(
  me = "#64a776",
  lp = "#537eb7",
  gc = "#d69971",
  vs = "#9f8d89",
  ep = "#b85292",
  fb = "#efd2c9"
)

# Step 21. Define full cell type names
cell_names <- c(
  me = "Mesophyll",
  lp = "Large parenchyma",
  gc = "Guard cell",
  vs = "Vascular",
  ep = "Epidermis",
  fb = "Fiber"
)

# Step 22. Plot correlations for all cell types
for (ct in celltypes) {

  df <- merged_list[[ct]] %>%
    filter(is.finite(rna_logFC), is.finite(peak_logFC))

  if (nrow(df) < 2 ||
      n_distinct(df$rna_logFC) < 2 ||
      n_distinct(df$peak_logFC) < 2) {
    warning(paste("Insufficient data for:", ct))
    next
  }

  cor_val <- cor(
    df$peak_logFC,
    df$rna_logFC,
    method = "spearman"
  )

  p <- ggplot(
    df,
    aes(x = peak_logFC, y = rna_logFC)
  ) +
    geom_point(
      size = 2,
      alpha = 0.7,
      color = cell_colors[[ct]]
    ) +
    geom_smooth(
      method = "lm",
      se = TRUE,
      color = "black",
      fill = "black",
      alpha = 0.25
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "grey70"
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey70"
    ) +
    theme_classic(base_size = 14) +
    labs(
      title = paste0(
        cell_names[[ct]],
        " (Spearman rho = ",
        round(cor_val, 2),
        ")"
      ),
      x = "Differential H3K27ac (avg_log2FC)",
      y = "RNA DEGs (avg_log2FC)"
    )

  ggsave(
    filename = file.path(
      output_dir,
      paste0("H3K27ac_RNA_correlation_", ct, ".pdf")
    ),
    plot = p,
    device = "pdf",
    width = 7,
    height = 6
  )

  print(p)
}
