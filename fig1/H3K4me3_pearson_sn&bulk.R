df <- read.table("E:/sn_CT/table/H3K4me3_readCounts.tab", header = FALSE)
colnames(df) <- c("chr", "start", "end", "H3K4me3_Nano_Tn5_rep1", "H3K4me3_Nano_Tn5_rep2", "H3K4me3_Nano_Tn5_Beads_rep1", "H3K4me3_Nano_Tn5_Beads_rep2", "H3K4me3_yl_rep1", "H3K4me3_yl_rep2", "H3K4me3_ml_rep1", "H3K4me3_ml_rep2")
signal <- df[, 4:ncol(df)]
cor_matrix <- cor(signal, method = "pearson")
library(pheatmap)

pdf("E:/sn_CT/figure/H3K4me3_correlation_heatmap.pdf", width = 8, height = 8)
pheatmap(
  cor_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,          # 显示相关系数数值
  number_color = "black",          # 数值字体颜色
  fontsize_number = 10,
  color = colorRampPalette(c("white", "#fdbb84", "#e34a33"))(100),  # 柔和、接近的颜色
  breaks = seq(0.85, 1, length.out = 101),  # 设置颜色区间，增强对比度
  border_color = NA,
  treeheight_row = 20,
  treeheight_col = 20
)
dev.off()



cor?
