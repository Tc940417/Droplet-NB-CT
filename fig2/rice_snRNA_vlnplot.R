color <- c('#927c9a','#3674a2','#9f8d89','#72567a',
           '#63a3b8','#c4daec','#61bada','#b7deea','#e29eaf',
           '#4490c4','#e6e2a3','#de8b36','#c4612f','#9a70a8',
           '#76a2be','#408444','#c6adb0','#9d3b62','#2d3462')
nip_filtered <- subset(nip, subset = cell_type != "Unknown")
desired_order <- c("Mesophyll", "Large parenchyma", "Guard cell", "Fiber", "Epidermis", "Vascular")
nip_filtered$cell_type <- factor(nip_filtered$cell_type, levels = desired_order)

# 你已有的基因列表
genes_to_check = c(
                   "LOC-Os05g02420",
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
                   "LOC-Os07g38960"
                   )

# 绘图
p1 <- VlnPlot(nip_filtered,
              features = genes_to_check,
              group.by = "cell_type",
              flip = TRUE,
              stack = TRUE, cols = color) + NoLegend() +
  theme(
    axis.text.x = element_text(size = 10),  # x轴字体大小
    axis.text.y = element_text(size = 3)    # y轴字体大小
  )
p1






