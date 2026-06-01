
seurat_list <- list(sample1_rep1 = ctrl_H3K4me3_Nano_Tn5_rep1,sample1_rep2 = ctrl_H3K4me3_Nano_Tn5_rep2,
                    sample2_rep1 = ctrl_H3K4me3_Nano_Tn5_Beads_rep1,sample2_rep2 = ctrl_H3K4me3_Nano_Tn5_Beads_rep2)
##TSS##

TSS_enrichment_df <- imap_dfr(seurat_list, ~{
  data.frame(
    cell_id = colnames(.x),
    TSS_enrichment = .x$TSS.enrichment,
    sample = .y,
    stringsAsFactors = FALSE
  )
}) %>% 
  tidyr::pivot_wider(
    names_from = sample,
    values_from = TSS_enrichment
  )


TSS_enrichment_df_1 <- TSS_enrichment_df[,-1]

colnames(TSS_enrichment_df_1) <- c("ctrl_H3K4me3_Nano_Tn5_rep1", "ctrl_H3K4me3_Nano_Tn5_rep2", 
                                    "ctrl_H3K4me3_Nano_Tn5_Beads_rep1","ctrl_H3K4me3_Nano_Tn5_Beads_rep2"
)

sample_colors <- c(
  rep('#df5734', 2),  # 第1-2列
  rep('#64a776', 2)),  # 第3-4列
  rep('#2d3462', 2),  # 第5-6列
  rep('#ebb1a4', 2),  # 第7-8列
  rep('#76a2be', 2)   # 第9-10列
)

plot_long <- TSS_enrichment_df_1 %>%
  pivot_longer(
    cols = everything(),
    names_to = "Sample",
    values_to = "TSS_enrichment"
  ) %>%
  mutate(
    Sample = factor(Sample, levels = colnames(TSS_enrichment_df_1)), # 保持原始列顺序
    Group = rep(1:2, each = 2)[as.numeric(Sample)] # 创建分组变量用于颜色映射
  )

plot <- ggplot(plot_long, aes(x = Sample, y = TSS_enrichment , fill = factor(Group))) +
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
    title = "TSS enrichment (Rice)",
    x = NULL,
    y = "TSS enrichment per cell",
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
    yintercept = c(1, 2), 
    linetype = "dashed", 
    color = "gray50",
    alpha = 0.5
  )
plot
ggsave(
  filename = "rice_broad_TSS_violin_plot.pdf",
  plot = plot,
  device = "pdf",
  width = 10,        # 宽度（英寸）
  height = 6,        # 高度（英寸）
  units = "in",      # 尺寸单位
  dpi = 1200          # 分辨率
)