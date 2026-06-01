sn_H3K4me3 <- readRDS("snCT/sn_H3K4me3.anno.rds")
sn_H3K27ac <- readRDS("snCT/sn_h3k27ac.anno.rds")

###gene.activities##
gene.activities1 <- GeneActivity(sn_H3K4me3)
gene.activities2 <- GeneActivity(sn_H3K27ac)

sn_H3K4me3[['RNA']] <- CreateAssayObject(counts = gene.activities1)
sn_H3K4me3 <- NormalizeData(
  object = sn_H3K4me3,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(sn_H3K4me3$nCount_RNA)
)

sn_H3K27ac[['RNA']] <- CreateAssayObject(counts = gene.activities2)
sn_H3K27ac <- NormalizeData(
  object = sn_H3K27ac,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(sn_H3K27ac$nCount_RNA)
)

DefaultAssay(sn_H3K4me3) <- "RNA"
DefaultAssay(sn_H3K27ac) <- "RNA"

Idents(sn_H3K4me3) <- "cell_type"
Idents(sn_H3K27ac) <- "cell_type"

saveRDS(sn_H3K4me3, file = "sn_H3K4me3.anno.addRNA.rds")
saveRDS(sn_H3K27ac, file = "sn_H3K27ac.anno.addRNA.rds")

sn_H3K4me3 <- readRDS("E:/sn_CT/RDS&h5ad/snCTaddRNA/sn_H3K4me3.anno.addRNA.rds")
sn_H3K27ac <- readRDS("E:/sn_CT/RDS&h5ad/snCTaddRNA/sn_H3K27ac.anno.addRNA.rds")
DefaultAssay(sn_H3K4me3) <- "RNA"
DefaultAssay(sn_H3K27ac) <- "RNA"

lable <- c("LOC-Os05g02420",
           "LOC-Os06g47640",
           "LOC-Os04g58800",
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
           "LOC-Os07g38960")

###sn_H3K27ac_avg-plot####
sn_H3K27ac_markers <- FindAllMarkers(
  sn_H3K27ac,
  assay = "RNA",
  only.pos = TRUE,         # 只考虑上调
  min.pct = 0.15,           # 表达比例阈值
  logfc.threshold = 0.15   # 初步阈值，可放宽
)

sn_H3K27ac_unique_markers <- sn_H3K27ac_markers %>%
  group_by(gene) %>%
  slice_max(order_by = avg_log2FC, n = 1) %>%
  ungroup()%>%
  arrange(cluster)


plot_order <- c(
  "Mesophyll",
  "Large parenchyma",
  "Epidermis",
  "Vascular",
  "Fiber",
  "Guard cell"
)
plot_order <- rev(plot_order)

sn_H3K27ac_unique_markers$cluster <- factor(
  sn_H3K27ac_unique_markers$cluster,
  levels = plot_order
)

sn_H3K27ac_unique_markers <- sn_H3K27ac_unique_markers[
  order(sn_H3K27ac_unique_markers$cluster),
]



sn_H3K27ac$cell_type_plot <- factor(
  sn_H3K27ac$cell_type,
  levels = plot_order
)

sn_H3K27ac <- sn_H3K27ac[, order(sn_H3K27ac$cell_type_plot)]





sn_H3K27ac_ht <- GroupHeatmap(
  srt = sn_H3K27ac,
  group.by = "cell_type_plot",
  features = sn_H3K27ac_unique_markers$gene,
  features_label = lable,
  height = 6.5,
  width = 4,
  group_palette = "simspec",
  heatmap_palcolor = c("#64a776", "white", "#537eb7"),
  column_title = "sn_H3K27ac Average expression heatmap"
)
sn_H3K27ac_ht

write.csv(
  sn_H3K27ac_markers,
  file = "E:/sn_CT/table/sn_H3K27ac_markers.csv",
  row.names = TRUE,
  quote = FALSE
)
###sn_H3K4me3_avg-plot####
sn_H3K4me3_markers <- FindAllMarkers(
  sn_H3K4me3,
  assay = "RNA",
  only.pos = TRUE,         # 只考虑上调
  min.pct = 0.25,           # 表达比例阈值
  logfc.threshold = 0.25   # 初步阈值，可放宽
)
sn_H3K4me3_unique_markers <- sn_H3K4me3_markers %>%
  group_by(gene) %>%
  slice_max(order_by = avg_log2FC, n = 1) %>%
  ungroup()%>%
  arrange(cluster)

sn_H3K4me3_ht <- GroupHeatmap(
  srt = sn_H3K4me3,
  group.by = "cell_type",
  features = sn_H3K4me3_unique_markers$gene,
  features_label = lable,
  height = 6.5,
  width = 4,
  group_palette = "simspec",
  heatmap_palcolor = c("#64a776", "white", "#537eb7"),
  column_title = "sn_H3K4me3 Average expression heatmap"
)
sn_H3K4me3_ht
write.csv(
  sn_H3K4me3_markers,
  file = "E:/sn_CT/table/sn_H3K4me3_markers.csv",
  row.names = TRUE,
  quote = FALSE
)
