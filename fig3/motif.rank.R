library(ggplot2)
library(dplyr)
df <- vs.da %>%
  tibble::rownames_to_column("motif") %>%
  arrange(desc(avg_diff)) %>%
  mutate(
    rank = row_number(),
    label = ifelse(rank <= 5 | rank > max(rank) - 5, motif, NA)
  )


ggplot(df, aes(x = rank, y = avg_diff, color = avg_diff)) +
  geom_point(size = 1.5) +
  geom_text_repel(
    aes(label = label),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey50"
  ) +
  scale_color_gradientn(
    colors = c("#053061", "#2166AC", "#F7F7F7", "#B2182B"),
    values = scales::rescale(c(min(df$avg_diff), -1, 0, max(df$avg_diff)))
  ) +
  theme_classic() +
  labs(
    x = "Rank (sorted by avg_diff)",
    y = "avg_diff (chromVAR deviation)",
    color = "avg_diff"
  )

