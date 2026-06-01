order <- c("mosophyll","phloem","myrosin_idioblast","late_S_phase/G2/M_phase","G2/M_phase",
           "late_S_phase","vascular","epidermis","guard_cell","companion_cell","xylem")


p <- DotPlot(Ara, assay = "RNA",
             features = c("AT3G63140", "AT2G05100", "AT4G13770", "AT5G25980",
                          "AT5G26000", "AT1G18370", "AT4G23800", "AT4G11080",
                          "AT5G54640", "AT4G20270", "AT5G43760", "AT1G25450",
                          "AT2G21140", "AT5G57350", "AT1G79430", "AT5G19530",
                          "AT4G32880", "AT5G60490"),
             group.by = "cell_type", cols = c("#ffffff", "firebrick3")) +
  scale_y_discrete(limits = rev(order)) +  # 👈 反转 y 轴顺序
  coord_flip() +
  theme(
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1, face = "italic"),
    axis.text.y = element_text(size = 8)
  ) +
  labs(x = NULL, y = NULL) +
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))

p
