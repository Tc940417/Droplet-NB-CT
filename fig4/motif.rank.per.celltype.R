
# ============================================================
# H3K4me3 Differential Motif Rank Comparison in Rice
# Spearman Rank Correlation and Alluvial Visualization
# ============================================================

# Step 1. Load required libraries
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(ggalluvial)

# Step 2. Set input and output directories
data_dir <- "results/H3K4me3_differential_motifs"
output_dir <- "results/H3K4me3_motif_rank"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Step 3. Define cell types and input files
motif_files <- c(
  Mesophyll = "H3K4me3_heat_vs_ctrl_Mesophyll_differential_motif_activity.csv",
  Large_parenchyma = "H3K4me3_heat_vs_ctrl_Large_parenchyma_differential_motif_activity.csv",
  Epidermis = "H3K4me3_heat_vs_ctrl_Epidermis_differential_motif_activity.csv",
  Vascular = "H3K4me3_heat_vs_ctrl_Vascular_differential_motif_activity.csv",
  Fiber = "H3K4me3_heat_vs_ctrl_Fiber_differential_motif_activity.csv",
  Guard_cell = "H3K4me3_heat_vs_ctrl_Guard_cell_differential_motif_activity.csv"
)

# Step 4. Read differential motif activity results
da_list <- lapply(motif_files, function(file) {
  read.csv(
    file.path(data_dir, file),
    row.names = 1,
    check.names = FALSE
  )
})

# Remove motifs with missing or invalid deviation differences
da_list <- lapply(da_list, function(df) {
  df[is.finite(df$avg_diff), , drop = FALSE]
})

# ============================================================
# Part 1. Identify Common Motifs Across Cell Types
# ============================================================

# Step 5. Identify motifs shared by all six cell types
common_motifs <- Reduce(
  intersect,
  lapply(da_list, rownames)
)

message("Number of common motifs: ", length(common_motifs))

if (length(common_motifs) < 2) {
  stop("Insufficient common motifs for rank correlation analysis.")
}

# Step 6. Rank common motifs within each cell type
rank_common_df <- lapply(names(da_list), function(ct) {

  da_list[[ct]] %>%
    rownames_to_column("feature") %>%
    filter(feature %in% common_motifs) %>%
    arrange(desc(avg_diff), feature) %>%
    mutate(
      cell_type = ct,
      rank = row_number()
    ) %>%
    select(feature, cell_type, rank, avg_diff)

}) %>%
  bind_rows()

# Step 7. Save common motif rankings
write.csv(
  rank_common_df,
  file = file.path(output_dir, "H3K4me3_common_motif_rankings.csv"),
  row.names = FALSE
)

# ============================================================
# Part 2. Spearman Rank Correlation Analysis
# ============================================================

# Step 8. Convert motif rankings to wide format
rank_common_wide <- rank_common_df %>%
  select(feature, cell_type, rank) %>%
  pivot_wider(
    names_from = cell_type,
    values_from = rank
  )

# Step 9. Create motif ranking matrix
rank_mat <- rank_common_wide %>%
  column_to_rownames("feature") %>%
  as.matrix()

# Step 10. Calculate Spearman rank correlations
cor_mat <- cor(
  rank_mat,
  method = "spearman",
  use = "pairwise.complete.obs"
)

# Step 11. Save rank correlations
write.csv(
  cor_mat,
  file = file.path(output_dir, "H3K4me3_motif_rank_correlation.csv"),
  row.names = TRUE,
  quote = FALSE
)

# ============================================================
# Part 3. Prepare Motif Rank Categories
# ============================================================

# Step 12. Define motif rank categories
rank_alluvial <- rank_common_df %>%
  mutate(
    rank_bin = case_when(
      rank <= 50  ~ "Top 50",
      rank <= 100 ~ "Top 100",
      rank <= 200 ~ "Top 200",
      rank <= 300 ~ "Top 300",
      rank <= 400 ~ "Top 400",
      TRUE        ~ ">400"
    )
  )

# Step 13. Set rank and cell type order
rank_alluvial$rank_bin <- factor(
  rank_alluvial$rank_bin,
  levels = c(
    "Top 50",
    "Top 100",
    "Top 200",
    "Top 300",
    "Top 400",
    ">400"
  )
)

rank_alluvial$cell_type <- factor(
  rank_alluvial$cell_type,
  levels = names(motif_files)
)

# Step 14. Save alluvial plot data
write.csv(
  rank_alluvial,
  file = file.path(output_dir, "H3K4me3_motif_rank_alluvial_data.csv"),
  row.names = FALSE
)

# ============================================================
# Part 4. Generate Motif Rank Alluvial Plot
# ============================================================

# Step 15. Define rank category colors
rank_colors <- c(
  "Top 50"  = "#76a2be",
  "Top 100" = "#ebb1a4",
  "Top 200" = "#2d3462",
  "Top 300" = "#ece399",
  "Top 400" = "#c376a7",
  ">400"    = "#64a776"
)

# Step 16. Plot motif rank correspondence
p <- ggplot(
  rank_alluvial,
  aes(
    x = cell_type,
    y = 1,
    stratum = rank_bin,
    alluvium = feature,
    fill = rank_bin
  )
) +
  geom_flow(
    stat = "alluvium",
    alpha = 0.3
  ) +
  geom_stratum(
    width = 0.4,
    color = "black"
  ) +
  scale_fill_manual(
    values = rank_colors
  ) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Number of motifs",
    fill = "Motif rank",
    title = "Rank-layer correspondence of H3K4me3 motifs across cell types"
  ) +
  theme(
    axis.text.x = element_text(
      size = 10,
      color = "black"
    ),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.position = "right"
  )

# Step 17. Save alluvial plot
ggsave(
  filename = file.path(
    output_dir,
    "H3K4me3_motif_rank_alluvial.pdf"
  ),
  plot = p,
  device = "pdf",
  width = 12,
  height = 7
)

# Step 18. Display results
print(cor_mat)
print(p)
