
# ============================================================
# Differential Motif Activity Analysis in Rice
# H3K4me3 chromVAR Motif Activity Across Cell Types
# ============================================================

# Step 1. Load required libraries
library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(ggrepel)

# Step 2. Set input and output directories
data_dir <- "data"
output_dir <- "results/chromVAR"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Step 3. Read H3K4me3 Seurat object
h3k4 <- readRDS(
  file.path(data_dir, "ct.h3k4.rds")
)

# Step 4. Set chromVAR assay and cell type identities
DefaultAssay(h3k4) <- "chromvar"
Idents(h3k4) <- "cell_type"

# Step 5. Identify differential motif activity across cell types
da.motif <- FindAllMarkers(
  object = h3k4,
  assay = "chromvar",
  test.use = "wilcox",
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  logfc.threshold = 0,
  min.pct = 0,
  only.pos = FALSE
)

# Step 6. Save differential motif activity results
write.csv(
  da.motif,
  file = file.path(output_dir, "H3K4me3_celltype_differential_motifs.csv"),
  row.names = FALSE
)

# Step 7. Rank motifs within each cell type
df <- da.motif %>%
  filter(is.finite(avg_diff)) %>%
  group_by(cluster) %>%
  arrange(desc(avg_diff), .by_group = TRUE) %>%
  mutate(
    rank = row_number(),
    label = ifelse(
      rank <= 5 | rank > n() - 5,
      gene,
      NA_character_
    )
  ) %>%
  ungroup()

# Step 8. Define motif activity color scale
motif_colors <- c(
  "#053061",
  "#2166AC",
  "#F7F7F7",
  "#B2182B"
)

# Step 9. Plot motif ranking for each cell type
for (ct in unique(df$cluster)) {

  plot_df <- df %>%
    filter(cluster == ct)

  p <- ggplot(
    plot_df,
    aes(x = rank, y = avg_diff, color = avg_diff)
  ) +
    geom_point(size = 1.5) +
    geom_text_repel(
      data = subset(plot_df, !is.na(label)),
      aes(label = label),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.4,
      point.padding = 0.3,
      segment.color = "grey50"
    ) +
    scale_color_gradientn(
      colours = motif_colors,
      values = scales::rescale(
        c(
          min(plot_df$avg_diff),
          -1,
          0,
          max(plot_df$avg_diff)
        )
      )
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "grey50"
    ) +
    theme_classic() +
    labs(
      title = paste("Differential Motif Activity -", ct),
      x = "Rank (sorted by avg_diff)",
      y = "avg_diff (chromVAR deviation)",
      color = "avg_diff"
    )

  # Save one plot per cell type
  ggsave(
    filename = file.path(
      output_dir,
      paste0(
        "H3K4me3_motif_ranking_",
        gsub("[^A-Za-z0-9]+", "_", ct),
        ".pdf"
      )
    ),
    plot = p,
    device = "pdf",
    width = 10,
    height = 6
  )

  print(p)
}

# Step 10. Save motif ranking table
write.csv(
  df,
  file = file.path(output_dir, "H3K4me3_celltype_motif_rankings.csv"),
  row.names = FALSE
)
