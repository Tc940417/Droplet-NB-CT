
# ============================================================
# Unique Fragments Analysis
# Comparison of Nano-Tn5 and Nano-Tn5-Beads in Rice
# ============================================================

# Step 1. Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)

# Step 2. Read Seurat objects
data_dir <- "data"
output_dir <- "results"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

seurat_list <- list(
  ctrl_H3K4me3_Nano_Tn5_rep1 = readRDS(file.path(data_dir, "snCT_5.rds")),
  ctrl_H3K4me3_Nano_Tn5_rep2 = readRDS(file.path(data_dir, "snCT_6.rds")),
  ctrl_H3K4me3_Nano_Tn5_Beads_rep1 = readRDS(file.path(data_dir, "snCT_7.rds")),
  ctrl_H3K4me3_Nano_Tn5_Beads_rep2 = readRDS(file.path(data_dir, "snCT_8.rds"))
)

# Step 3. Extract log10 unique fragments
frags_df <- bind_rows(lapply(names(seurat_list), function(sample) {
  data.frame(
    Sample = sample,
    log10_uniqueFrags = seurat_list[[sample]]$log10_uniqueFrags
  )
}))

# Remove missing values
frags_df <- frags_df %>%
  filter(!is.na(log10_uniqueFrags))

# Step 4. Define sample groups and plotting order
frags_df$Sample <- factor(
  frags_df$Sample,
  levels = names(seurat_list)
)

frags_df$Group <- ifelse(
  grepl("Beads", frags_df$Sample),
  "Nano-Tn5-Beads",
  "Nano-Tn5"
)

frags_df$Group <- factor(
  frags_df$Group,
  levels = c("Nano-Tn5", "Nano-Tn5-Beads")
)

# Step 5. Plot unique fragments distributions
p <- ggplot(frags_df, aes(x = Sample, y = log10_uniqueFrags, fill = Group)) +
  geom_violin(
    width = 0.8,
    alpha = 0.7,
    scale = "width",
    trim = TRUE
  ) +
  geom_boxplot(
    width = 0.15,
    fill = "white",
    outlier.shape = NA
  ) +
  scale_fill_manual(
    values = c(
      "Nano-Tn5" = "#df5734",
      "Nano-Tn5-Beads" = "#64a776"
    ),
    labels = c(
      "ctrl_H3K4me3_Nano_Tn5",
      "ctrl_H3K4me3_Nano_Tn5_Beads"
    )
  ) +
  geom_hline(
    yintercept = c(3.5, 4.5, 5.5),
    linetype = "dashed",
    color = "gray50",
    alpha = 0.5
  ) +
  labs(
    title = "Unique Fragments (Rice)",
    x = NULL,
    y = expression("Unique Fragments per cell (log"[10]*")"),
    fill = "Sample Group"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# Step 6. Save the figure
ggsave(
  filename = file.path(output_dir, "rice_broad_Unique_Fragments_violin_plot.pdf"),
  plot = p,
  device = "pdf",
  width = 10,
  height = 6,
  units = "in"
)

print(p)
