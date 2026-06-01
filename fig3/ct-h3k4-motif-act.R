h3k4 <- readRDS("E:/sn_CT/RDS&h5ad/snCTaddRNAaddChromvar/ct.h3k4.rds")
motif_avg <- AverageExpression(
  h3k4,
  assays = "chromvar",
  slot = "data",
  group.by = "cell_type"
)
motif_avg_matrix <- motif_avg$chromvar
motif_avg_matrix <- as.matrix(motif_avg$chromvar)

library(ComplexHeatmap)
library(circlize)

Heatmap(
  motif_avg_matrix,
  name = "Deviation",
  col = colorRamp2(
    c(-2, 0, 2),
    c("#f7f7f7", "#c7e9e6", "#2b5c8a")
  ),
  show_row_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE
)


top10_list <- lapply(
  colnames(motif_avg_matrix),
  function(ct) {
    vals <- motif_avg_matrix[, ct]
    top_motifs <- names(sort(vals, decreasing = TRUE))[1:10]
    data.frame(
      cell_type = ct,
      motif = top_motifs,
      deviation = vals[top_motifs]
    )
  }
)

top10_motif_df <- do.call(rbind, top10_list)
top10_motif_df

df_split <- split(top10_motif_df, top10_motif_df$cell_type)

MotifPlot(
  object = h3k4,
  motifs = df_split$Vascular$motif,
  assay = 'ATAC'
)


library(dplyr)
library(tidyr)

df_long <- motif_avg_matrix %>%
  as.data.frame() %>%
  tibble::rownames_to_column("motif") %>%
  pivot_longer(
    cols = -motif,
    names_to = "cell_type",
    values_to = "avg_dev"
  )

df_long <- df_long %>%
  group_by(cell_type) %>%
  arrange(desc(avg_dev)) %>%
  mutate(rank = row_number()) %>%
  ungroup()
df_long <- df_long %>%
  group_by(cell_type) %>%
  mutate(
    rank_scaled = (rank - min(rank)) / (max(rank) - min(rank))
  ) %>%
  ungroup()
df_long <- df_long %>%
  group_by(cell_type) %>%
  mutate(
    label = ifelse(rank <= 5, motif, NA)
  ) %>%
  ungroup()

library(ggplot2)
library(ggrepel)

ggplot(
  df_long %>% filter(cell_type == "Fiber"),
  aes(x = rank, y = avg_dev, color = rank_scaled)
) +
  geom_point(size = 1.4) +
  geom_text_repel(
    aes(label = label),
    nudge_x = 320,
    direction = "y",
    hjust = 0,
    size = 3,
    fontface = "bold",
    max.overlaps = Inf,
    segment.color = "grey50",
    segment.curvature = 0
  ) +
  scale_color_gradientn(
    colors = c("#D7191C", "#FDB863", "#2C7BB6"),
    values = c(0, 0.5, 1),
    name = "TF rank",
    guide = guide_colorbar(reverse = TRUE)  # 👈 关键
  ) +
  scale_x_continuous(
    limits = c(min(df_long$rank), max(df_long$rank) + 0)
  ) +
  theme_classic() +
  labs(
    x = "TF rank",
    y = "Mean chromVAR deviation"
  )




write.csv(
  motif_avg_matrix,
  file = "E:/sn_CT/table/h3k4motif_avg_matrix.csv",
  row.names = TRUE,
  quote = FALSE
)









