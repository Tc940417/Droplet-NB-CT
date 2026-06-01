order <- c("Mesophyll","Large parenchyma","Guard cell","Fiber","Epidermis","Vascular")
p <- DotPlot(nip, assay = "RNA",
             features = c("LOC-Os05g02420",
                          "LOC-Os06g47640",
                          "LOC-Os04g58800",
                          "LOC-Os03g55240",
                          "LOC-Os02g41860",
                          "LOC-Os04g20164",
                          "LOC-Os12g02320",
                          "LOC-Os11g10310",
                          "LOC-Os11g10320",
                          "LOC-Os09g10054",
                          "LOC-Os08g43670",
                          "LOC-Os07g46852",
                          "LOC-Os01g15010",
                          "LOC-Os10g42190",
                          "LOC-Os01g64730",
                          "LOC-Os11g13890",
                          "LOC-Os09g17740",
                          "LOC-Os08g33820",
                          "LOC-Os07g38960"),
             group.by = "cell_type", cols = c("#ffffff", "firebrick3")) +
  scale_y_discrete(limits = factor(order)) +
  coord_flip() +
  theme(
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1, face = "italic"),
    axis.text.y = element_text(size = 8)  # 👈 调整 y 轴字体大小
  ) +
  labs(x = NULL, y = NULL) +
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))
p
