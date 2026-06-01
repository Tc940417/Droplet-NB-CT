ctrl_H3K4me3_Nano_Tn5_rep1 <- snCT_5
ctrl_H3K4me3_Nano_Tn5_rep2 <- snCT_6
ctrl_H3K4me3_Nano_Tn5_Beads_rep1 <- snCT_7
ctrl_H3K4me3_Nano_Tn5_Beads_rep2 <- snCT_8
ctrl_H3K27ac_Nano_Tn5_Beads_rep1 <- snCT_13
ctrl_H3K27ac_Nano_Tn5_Beads_rep2 <- snCT_14
heat_H3K4me3_Nano_Tn5_Beads_rep1 <- snCT_15
heat_H3K4me3_Nano_Tn5_Beads_rep2 <- snCT_16
heat_H3K27ac_Nano_Tn5_Beads_rep1 <- snCT_17
heat_H3K27ac_Nano_Tn5_Beads_rep2 <- snCT_18





seurat_list <- list(sample1_rep1 = ctrl_H3K4me3_Nano_Tn5_rep1,sample1_rep2 = ctrl_H3K4me3_Nano_Tn5_rep2,
                    sample2_rep1 = ctrl_H3K4me3_Nano_Tn5_Beads_rep1,sample2_rep2 = ctrl_H3K4me3_Nano_Tn5_Beads_rep2),
                    sample3_rep1 = ctrl_H3K27ac_Nano_Tn5_Beads_rep1,sample3_rep2 = ctrl_H3K27ac_Nano_Tn5_Beads_rep2,
                    sample4_rep1 = heat_H3K4me3_Nano_Tn5_Beads_rep1,sample4_rep2 = heat_H3K4me3_Nano_Tn5_Beads_rep2,
                    sample5_rep1 = heat_H3K27ac_Nano_Tn5_Beads_rep1,sample5_rep2 = heat_H3K27ac_Nano_Tn5_Beads_rep2)

##FRiP##
FRiP_df <- imap_dfr(seurat_list, ~{
  data.frame(
    cell_id = colnames(.x),
    FRiP = .x$FRiP,
    sample = .y,
    stringsAsFactors = FALSE
  )
}) %>% 
  tidyr::pivot_wider(
    names_from = sample,
    values_from = FRiP
  )


frip_df_1 <- FRiP_df[,-1]

colnames(frip_df_1) <- c("ctrl_H3K4me3_Nano_Tn5_rep1", "ctrl_H3K4me3_Nano_Tn5_rep2", 
                         "ctrl_H3K4me3_Nano_Tn5_Beads_rep1","ctrl_H3K4me3_Nano_Tn5_Beads_rep2"), 
#                         "ctrl_H3K27ac_Nano_Tn5_Beads_rep1", "ctrl_H3K27ac_Nano_Tn5_Beads_rep2", 
#                         "heat_H3K4me3_Nano_Tn5_Beads_rep1", "heat_H3K4me3_Nano_Tn5_Beads_rep2", 
#                         "heat_H3K27ac_Nano_Tn5_Beads_rep1", "heat_H3K27ac_Nano_Tn5_Beads_rep2" 
#)

sample_colors <- c(
  rep('#df5734', 2),  # 第1-2列
  rep('#64a776', 2)),  # 第3-4列
  rep('#2d3462', 2),  # 第5-6列
  rep('#ebb1a4', 2),  # 第7-8列
  rep('#76a2be', 2)   # 第9-10列
)

frip_long <- frip_df_1 %>%
  pivot_longer(
    cols = everything(),
    names_to = "Sample",
    values_to = "FRiP"
  ) %>%
  mutate(
    Sample = factor(Sample, levels = colnames(frip_df_1)), # 保持原始列顺序
    Group = rep(1:2, each = 2)[as.numeric(Sample)] # 创建分组变量用于颜色映射
  )

plot <- ggplot(frip_long, aes(x = Sample, y = FRiP, fill = factor(Group))) +
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
    values = c('#df5734','#64a776',
               '#2d3462','#ebb1a4','#76a2be'),
    labels = c(
      "ctrl_H3K4me3_Nano_Tn5",
      "ctrl_H3K4me3_Nano_Tn5_Beads",
      "ctrl_H3K27ac_Nano_Tn5_Beads",
      "heat_H3K4me3_Nano_Tn5_Beads",
      "heat_H3K27ac_Nano_Tn5_Beads"
    )
  ) +
  labs(
    title = "FRiP Distribution (Rice)",
    x = NULL,
    y = "Fraction of Reads in Peaks (FRiP)",
    fill = "Sample Group"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  geom_hline(
    yintercept = c(0.3, 0.6, 0.9), 
    linetype = "dashed", 
    color = "gray50",
    alpha = 0.5
  )
plot

ggsave(
  filename = "rice_broad_FRiP_violin_plot.pdf",
  plot = plot,
  device = "pdf",
  width = 10,        # 宽度（英寸）
  height = 6,        # 高度（英寸）
  units = "in",      # 尺寸单位
  dpi = 300          # 分辨率
)