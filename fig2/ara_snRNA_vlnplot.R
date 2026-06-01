color <- c('#4b6aa8','#3ca0cf','#c376a7','#ad98c3','#cea5c7',
           '#53738c','#a5a9b0','#a78982','#696a6c','#92699e',
           '#d69971','#df5734','#6c408e','#ac6894','#d4c2db',
           '#537eb7','#83ab8e','#ece399','#405993')

desired_order <- c("mosophyll","phloem","myrosin_idioblast","late_S_phase/G2/M_phase","G2/M_phase",
                   "late_S_phase","vascular","epidermis","guard_cell","companion_cell","xylem")
Araun$cell_type <- factor(Araun$cell_type, levels = desired_order)

# 你已有的基因列表
genes_to_check = c("AT3G63140", "AT2G05100", "AT4G13770", "AT5G25980",
                   "AT5G26000", "AT1G18370", "AT4G23800", "AT4G11080",
                   "AT5G54640", "AT4G20270", "AT5G43760", "AT1G25450",
                   "AT2G21140", "AT5G57350", "AT1G79430", "AT5G19530",
                   "AT4G32880", "AT5G60490")

genes_to_check_rev <- rev(genes_to_check)
# 绘图
p1 <- VlnPlot(Araun,
              features = genes_to_check_rev,
              group.by = "cell_type",
              flip = TRUE,
              stack = TRUE, cols = color) + NoLegend() +
  theme(
    axis.text.x = element_text(size = 10),  # x轴字体大小
    axis.text.y = element_text(size = 3)    # y轴字体大小
  )
p1
