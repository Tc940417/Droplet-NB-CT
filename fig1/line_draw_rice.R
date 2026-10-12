
# ============================================================
# Sequencing Saturation Analysis in Rice
# Mean Detected Fragments and Peaks vs. Sequencing Depth
# ============================================================

# Step 1. Load required libraries
library(ggplot2)
library(reshape2)

# Step 2. Set input and output directories
data_dir <- "Sequencing saturation"
output_dir <- "Sequencing saturation"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Step 3. Define sample colors
custom_colors <- c(
  "#df5734",
  "#e97e61",
  "#64a776",
  "#8fc19a"
)

# Step 4. Plot mean detected fragments vs. sequencing depth
count <- read.csv(
  file.path(data_dir, "means_nCount_by_nReads_custom.csv")
)

count_long <- melt(
  count,
  id.vars = "X",
  variable.name = "Group",
  value.name = "Value"
)

p1 <- ggplot(count_long, aes(x = X, y = Value, color = Group)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  scale_color_manual(values = custom_colors) +
  labs(
    title = "Mean of detected fragments",
    x = "Number of reads",
    y = "Mean of detected fragments",
    color = "Group"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 90, hjust = 1, color = "black",
      face = "bold", size = 12
    ),
    axis.text.y = element_text(
      color = "black", face = "bold", size = 12
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      color = "black", linewidth = 1.2, linetype = "solid"
    ),
    panel.spacing = grid::unit(0.12, "cm"),
    plot.title = element_text(
      hjust = 0.5, size = 14, face = "bold"
    )
  )

ggsave(
  filename = file.path(output_dir, "fragments-reads_rice.pdf"),
  plot = p1,
  device = "pdf",
  width = 8,
  height = 5,
  units = "in"
)

# Step 5. Plot mean detected peaks vs. sequencing depth
peak <- read.csv(
  file.path(data_dir, "means_nFeature_by_nReads_custom.csv")
)

peak_long <- melt(
  peak,
  id.vars = "X",
  variable.name = "Group",
  value.name = "Value"
)

p2 <- ggplot(peak_long, aes(x = X, y = Value, color = Group)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  scale_color_manual(values = custom_colors) +
  labs(
    title = "Mean of detected peaks",
    x = "Number of reads",
    y = "Mean of detected peaks",
    color = "Group"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 90, hjust = 1, color = "black",
      face = "bold", size = 12
    ),
    axis.text.y = element_text(
      color = "black", face = "bold", size = 12
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      color = "black", linewidth = 1.2, linetype = "solid"
    ),
    panel.spacing = grid::unit(0.12, "cm"),
    plot.title = element_text(
      hjust = 0.5, size = 14, face = "bold"
    )
  )

ggsave(
  filename = file.path(output_dir, "peaks-reads_rice.pdf"),
  plot = p2,
  device = "pdf",
  width = 8,
  height = 5,
  units = "in"
)

# Step 6. Display plots
print(p1)
print(p2)
