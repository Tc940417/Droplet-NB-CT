
# ============================================================
# chromVAR Motif Activity Analysis in Rice
# H3K4me3 Motif Heatmap, Top Motifs, and TF Ranking
# ============================================================

# Step 1. Load required libraries
library(Seurat)
library(Signac)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)
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

# ============================================================
# Part 1. Average chromVAR Motif Deviation Heatmap
# ============================================================

# Step 4. Extract chromVAR deviation scores
motif_matrix <- GetAssayData(
  h3k4,
  assay = "chromvar",
  layer = "data"
)

# Step 5. Calculate mean motif deviation per cell type
cell_types <- h3k4$cell_type
names(cell_types) <- colnames(h3k4)

valid_cells <- colnames(motif_matrix)[
  !is.na(cell_types[colnames(motif_matrix)])
]

motif_matrix <- motif_matrix[, valid_cells, drop = FALSE]
cell_types <- cell_types[valid_cells]

motif_avg_matrix <- sapply(
  unique(as.character(cell_types)),
  function(ct) {
    cells <- names(cell_types)[as.character(cell_types) == ct]
    Matrix::rowMeans(motif_matrix[, cells, drop = FALSE])
  }
)

motif_avg_matrix <- as.matrix(motif_avg_matrix)

# Step 6. Generate motif deviation heatmap
ht <- Heatmap(
  motif_avg_matrix,
  name = "Deviation",
  col = colorRamp2(
    c(-2, 0, 2),
    c("#f7f7f7", "#c7e9e6", "#2b5c8a")
  ),
  show_row_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_title = "H3K4me3 Motif Activity in Rice"
)

# Step 7. Save heatmap
pdf(
  file.path(output_dir, "H3K4me3_chromVAR_heatmap.pdf"),
  width = 8,
  height = 10
)

draw(ht)
dev.off()

# ============================================================
# Part 2. Identify Top 10 Motifs for Each Cell Type
# ============================================================

# Step 8. Extract top 10 motifs
top10_list <- lapply(
  colnames(motif_avg_matrix),
  function(ct) {

    vals <- motif_avg_matrix[, ct]
    vals <- sort(vals[is.finite(vals)], decreasing = TRUE)

    top_motifs <- head(names(vals), 10)

    data.frame(
      cell_type = ct,
      motif = top_motifs,
      deviation = unname(vals[top_motifs]),
      row.names = NULL
    )
  }
)

top10_motif_df <- bind_rows(top10_list)

# Step 9. Save average deviation matrix and top motifs
write.csv(
  motif_avg_matrix,
  file = file.path(output_dir, "H3K4me3_motif_avg_matrix.csv"),
  row.names = TRUE,
  quote = FALSE
)

write.csv(
  top10_motif_df,
  file = file.path(output_dir, "H3K4me3_top10_motifs.csv"),
  row.names = FALSE,
  quote = FALSE
)

# ============================================================
# Part 3. Plot Top Motif Sequence Logos in Vascular Cells
# ============================================================

# Step 10. Extract top motifs in Vascular cells
vascular_motifs <- top10_motif_df %>%
  filter(cell_type == "Vascular") %>%
  pull(motif)

# Step 11. Match motifs with the ATAC assay
available_motifs <- colnames(Motifs(h3k4[["ATAC"]]))

vascular_motifs <- intersect(
  vascular_motifs,
  available_motifs
)

# Step 12. Generate motif sequence logos
if (length(vascular_motifs) > 0) {

  p_logo <- MotifPlot(
    object = h3k4,
    motifs = vascular_motifs,
    assay = "ATAC"
  )

  ggsave(
    filename = file.path(
      output_dir,
      "H3K4me3_Vascular_top10_motif_logos.pdf"
    ),
    plot = p_logo,
    device = "pdf",
    width = 12,
    height = 8
  )

  print(p_logo)

} else {
  warning("No matching Vascular motifs found in the ATAC assay.")
}

# ============================================================
# Part 4. Rank Motif Activity in Fiber Cells
# ============================================================

# Step 13. Convert average deviation matrix to long format
df_long <- motif_avg_matrix %>%
  as.data.frame() %>%
  tibble::rownames_to_column("motif") %>%
  pivot_longer(
    cols = -motif,
    names_to = "cell_type",
    values_to = "avg_dev"
  )

# Step 14. Rank motifs within each cell type
df_long <- df_long %>%
  filter(is.finite(avg_dev)) %>%
  group_by(cell_type) %>%
  arrange(desc(avg_dev), motif, .by_group = TRUE) %>%
  mutate(
    rank = row_number(),
    rank_scaled = if (n() > 1) {
      (rank - 1) / (n() - 1)
    } else {
      0
    },
    label = ifelse(rank <= 5, motif, NA_character_)
  ) %>%
  ungroup()

# Step 15. Extract Fiber motif rankings
fiber_df <- df_long %>%
  filter(cell_type == "Fiber")

# Step 16. Plot Fiber motif rankings
p_rank <- ggplot(
  fiber_df,
  aes(x = rank, y = avg_dev, color = rank_scaled)
) +
  geom_point(size = 1.4) +
  geom_text_repel(
    aes(label = label),
    data = subset(fiber_df, rank <= 5),
    nudge_x = 320,
    direction = "y",
    hjust = 0,
    size = 3,
    fontface = "bold",
    max.overlaps = Inf,
    segment.color = "grey50"
  ) +
  scale_color_gradientn(
    colours = c("#D7191C", "#FDB863", "#2C7BB6"),
    values = c(0, 0.5, 1),
    name = "TF rank",
    guide = guide_colorbar(reverse = TRUE)
  ) +
  scale_x_continuous(
    limits = c(1, max(fiber_df$rank) + 350)
  ) +
  theme_classic() +
  labs(
    title = "H3K4me3 Motif Ranking in Fiber Cells",
    x = "TF rank",
    y = "Mean chromVAR deviation"
  )

# Step 17. Save Fiber motif ranking plot
ggsave(
  filename = file.path(
    output_dir,
    "H3K4me3_Fiber_motif_ranking.pdf"
  ),
  plot = p_rank,
  device = "pdf",
  width = 10,
  height = 6
)

# Step 18. Save complete motif rankings
write.csv(
  df_long,
  file = file.path(output_dir, "H3K4me3_motif_rankings.csv"),
  row.names = FALSE,
  quote = FALSE
)

# Step 19. Display plots
draw(ht)
print(p_rank)
