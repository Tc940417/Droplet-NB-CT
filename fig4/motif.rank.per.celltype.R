library(dplyr)
library(tibble)
library(tidyr)

da_list <- list(
  Mesophyll = me.da,
  Large_parenchyma = lp.da,
  Epidermis = ep.da,
  Vascular = vs.da,
  Fiber = fb.da,
  Guard_cell = gc.da
)
common_motifs <- Reduce(
  intersect,
  lapply(da_list, rownames)
)

length(common_motifs)
rank_common_df <- lapply(names(da_list), function(ct) {
  da_list[[ct]] %>%
    rownames_to_column("feature") %>%
    filter(feature %in% common_motifs) %>%
    arrange(desc(avg_diff)) %>%
    mutate(
      cell_type = ct,
      rank = row_number()
    ) %>%
    select(feature, cell_type, rank, avg_diff)
}) %>%
  bind_rows()

rank_common_wide <- rank_common_df %>%
  select(feature, cell_type, rank) %>%
  pivot_wider(
    names_from = cell_type,
    values_from = rank
  )

rank_mat <- rank_common_wide %>%
  column_to_rownames("feature") %>%
  as.matrix()

cor_mat <- cor(
  rank_mat,
  method = "spearman",
  use = "pairwise.complete.obs"
)

cor_mat




rank_alluvial <- rank_common_df %>%
  group_by(cell_type) %>%
  mutate(
    rank_bin = case_when(
      rank <= 50  ~ "Top 50",
      rank <= 100 ~ "Top 100",
      rank <= 200 ~ "Top 200",
      rank <= 300 ~ "Top 300",
      rank <= 400 ~ "Top 400",
      TRUE        ~ ">400"
    )
  ) %>%
  ungroup()
rank_alluvial <- rank_alluvial %>%
  mutate(
    rank_bin = factor(
      rank_bin,
      levels = c(
        "Top 50",
        "Top 100",
        "Top 200",
        "Top 300",
        "Top 400",
        ">400"
      )
    )
  )

library(ggalluvial)
library(ggplot2)

ggplot(
  rank_alluvial,
  aes(
    x = cell_type,
    stratum = rank_bin,
    alluvium = feature,
    fill = rank_bin
  )
) +
  geom_flow(
    stat = "alluvium",
    alpha = 0.3
  ) +
  geom_stratum(width = 0.4, color = "black") +
  scale_fill_manual(
    values = c(
      "Top 50"  = '#76a2be',
      "Top 100" = '#ebb1a4',
      "Top 200" = '#2d3462',
      "Top 300" = '#ece399',
      "Top 400" = '#c376a7',
      ">400"    = '#64a776'
    )
  ) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Number of motifs",
    title = "Rank-layer correspondence of H3K4me3 motifs across cell types"
  )