
# ============================================================
# Mitochondrial Percentage Analysis
# Comparison of Nano-Tn5 and Nano-Tn5-Beads
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

# Step 3. Extract mitochondrial percentages
mt_df <- bind_rows(lapply(names(seurat_list), function(sample) {
  data.frame(
    Sample = sample,
    percent.mt = seurat_list[[sample]]$percent.mt
  )
}))

# Step 4. Define sample groups and plotting order
mt_df$Sample <- factor(
  mt_df$Sample,
  levels = names(seurat_list)
)

mt_df$Group <- ifelse(
  grepl("Beads", mt_df$Sample),
  "Nano-Tn5-Beads",
  "Nano-Tn5"
)

mt_df$Group <- factor(
  mt_df$Group,
  levels = c("Nano-Tn5", "Nano-Tn5-Beads")
)

# Step 5. Plot mitochondrial percentages
p <- ggplot(mt_df, aes(x = Sample, y = percent.mt, fill = Group)) +
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
    yintercept = c(1, 3, 5),
    linetype = "dashed",
    color = "gray50",
    alpha = 0.5
  ) +
  labs(
    title = "Mitochondrial percentage",
    x = NULL,
    y = NULL,
    fill = "Sample Group"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# Step 6. Save plot
ggsave(
  filename = file.path(output_dir, "rice_Mt.percent_violin_plot.pdf"),
  plot = p,
  device = "pdf",
  width = 10,
  height = 6,
  units = "in"
)

print(p)
