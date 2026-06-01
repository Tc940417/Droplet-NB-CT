seurat_list <- list(sample2_rep1 = ctrl_H3K4me3_Nano_Tn5_rep1,sample2_rep2 = ctrl_H3K4me3_Nano_Tn5_rep2,
                    sample3_rep1 = ctrl_H3K4me3_Nano_Tn5_Beads_rep1,sample3_rep2 = ctrl_H3K4me3_Nano_Tn5_Beads_rep2,
                    sample4_rep1 = ctrl_H3K27ac_Nano_Tn5_Beads_rep1,sample4_rep2 = ctrl_H3K27ac_Nano_Tn5_Beads_rep2,
                    sample5_rep1 = heat_H3K4me3_Nano_Tn5_Beads_rep1,sample5_rep2 = heat_H3K4me3_Nano_Tn5_Beads_rep2,
                    sample6_rep1 = heat_H3K27ac_Nano_Tn5_Beads_rep1,sample6_rep2 = heat_H3K27ac_Nano_Tn5_Beads_rep2)




##unique_frags##

log10_uniqueFrags_df <- imap_dfr(seurat_list, ~{
  data.frame(
    cell_id = colnames(.x),
    log10_uniqueFrags = .x$log10_uniqueFrags,
    sample = .y,
    stringsAsFactors = FALSE
  )
}) %>% 
  tidyr::pivot_wider(
    names_from = sample,
    values_from = log10_uniqueFrags
  )


log10_uniqueFrags_df_1 <- log10_uniqueFrags_df[,-1]

colnames(log10_uniqueFrags_df_1) <- c("ctrl_H3K4me3_Nano_Tn5_rep1", "ctrl_H3K4me3_Nano_Tn5_rep2", 
                                      "ctrl_H3K4me3_Nano_Tn5_Beads_rep1","ctrl_H3K4me3_Nano_Tn5_Beads_rep2"
)

sample_colors <- c(
  rep('#df5734', 2),  # 
  rep('#64a776', 2)),  # 
  rep('#2d3462', 2),  # 
  rep('#ebb1a4', 2),  # 
  rep('#76a2be', 2)   # 
)

plot_long <- log10_uniqueFrags_df_1 %>%
  pivot_longer(
    cols = everything(),
    names_to = "Sample",
    values_to = "log10_uniqueFrags"
  ) %>%
  mutate(
    Sample = factor(Sample, levels = colnames(log10_uniqueFrags_df_1)), # 保持原始列顺序
    Group = rep(1:2, each = 2)[as.numeric(Sample)] # 创建分组变量用于颜色映射
  )

plot <- ggplot(plot_long, aes(x = Sample, y = log10_uniqueFrags , fill = factor(Group))) +
  # 小提琴图层
  geom_violin(
    width = 0.8, 
    alpha = 0.7,
    scale = "width",
    trim = TRUE
  ) +
  # 箱线图展示统计量
  geom_boxplot(
    width = 0.15,
    fill = "white",
    outlier.shape = NA
  ) +
  # 使用指定配色方案
  scale_fill_manual(
    values = c('#df5734','#64a776'),
    labels = c(
      "ctrl_H3K4me3_Nano_Tn5",
      "ctrl_H3K4me3_Nano_Tn5_Beads"
    )
  ) +
  # 美化主题
  labs(
    title = "Unique Fragments (Rice)",
    x = NULL,
    y = expression("Unique Fragments per cell (log"["10"]*")"),
    fill = "Sample Group"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  # 添加水平参考线
  geom_hline(
    yintercept = c(3.5, 4.5, 5.5), 
    linetype = "dashed", 
    color = "gray50",
    alpha = 0.5
  )
plot
ggsave(
  filename = "rice_broad_Unique_Fragments_violin_plot.pdf",
  plot = plot,
  device = "pdf",
  width = 10,        # 宽度（英寸）
  height = 6,        # 高度（英寸）
  units = "in",      # 尺寸单位
  dpi = 1200          # 分辨率
)