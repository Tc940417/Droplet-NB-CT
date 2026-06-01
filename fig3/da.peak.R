library(ggplot2)
library(tidyverse)
library(reshape2)
library(Seurat)
library(tidyr)
library(dplyr)
library(SCP)
library(BiocParallel)
library(ggplot2)
library(gtable)
library(Signac)
h3k4 <- readRDS("E:/sn_CT/RDS&h5ad/snCTaddRNAaddChromvar/ct.h3k4.rds")
DefaultAssay(h3k4) <- "ATAC"
h3k4 <- RunTFIDF(h3k4)
#h3k4.da <- FindAllMarkers(
#  h3k4,
#  test.use = "LR",
#  latent.vars = "nCount_ATAC",
#  min.pct = 0.1,
#  logfc.threshold = 0.25,
#  only.pos = TRUE
#)
#write.csv(
#  h3k4.da,
#  file = "E:/sn_CT/table/snCT_celltype_marker/ctH3K27_FindAllMarkers_LR.csv",
#  row.names = FALSE
#)

h3k4.da <- read.csv(
  "E:/sn_CT/table/snCT_celltype_marker/ctH3K4_FindAllMarkers_LR.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)





cf <- ClosestFeature(h3k4, regions = h3k4.da$gene)

h3k4.da <- cbind(h3k4.da, 
                       gene=cf$gene_name, 
                       query_region = cf$query_region, 
                       type = cf$type, 
                       distance=cf$distance)
h3k4.da <- h3k4.da[, -7]

#marker <- h3k4.da %>% dplyr::filter(avg_log2FC >=0.25 & p_val_adj < 0.05) %>% arrange(desc(avg_log2FC))
#sig.region <- marker %>% dplyr::select(query_region) %>% distinct()


sig.region.mean <- AverageExpression(h3k4,
                                     features = h3k4.da$query_region,
                                     assays = "ATAC")

sig.region.mean.scale <- scale(t(sig.region.mean$ATAC))




library(ComplexHeatmap)
library(circlize)
library(grid)

#==============================
# 1. 行顺序（细胞类型顺序）
#==============================
row_order <- c(
  "Epidermis",
  "Guard cell",
  "Mesophyll",
  "Large parenchyma",
  "Vascular",
  "Fiber"
)

sig.region.mean.scale <- sig.region.mean.scale[row_order, ]

#==============================
# 2. 列顺序（按主导细胞类型）
#==============================
col_order <- order(apply(sig.region.mean.scale, 2, which.max))
sig.region.mean.scale <- sig.region.mean.scale[, col_order]

# 每一列归属的细胞类型（用于标尺/分组）
col_group <- apply(sig.region.mean.scale, 2, function(x) {
  rownames(sig.region.mean.scale)[which.max(x)]
})

#==============================
# 3. 颜色（z-score）
#==============================
lim <- 2
col_fun <- circlize::colorRamp2(
  c(-lim, -1, 0, 1, lim),
  c("#6F8F7A",  # -2 灰绿（莫兰迪绿）
    "#B8CFC1",  # -1 浅灰绿
    "#FFFFFF",  #  0 白
    "#B7C9D9",  # +1 浅灰蓝
    "#5B7FA6")  # +2 灰蓝
)

#==============================
# 4. Heatmap
#==============================
ht <- Heatmap(
  sig.region.mean.scale,
  name = "z-score",
  col = col_fun,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  column_split = factor(col_group, levels = row_order),
  column_title = NULL,
  
  show_column_names = FALSE,
  
  row_names_gp = gpar(fontsize = 8, fontface = "bold"),
  
  width = unit(12, "cm"),
  height = unit(8, "cm"),
  
  # 🔹 黑色外边框
  border = TRUE,
  
  heatmap_legend_param = list(
    title = "Z-score",
    at = c(-2, -1, 0, 1, 2),
    labels = c("-2", "-1", "0", "1", "2"),
    title_gp = gpar(fontsize = 10),
    labels_gp = gpar(fontsize = 9),
    legend_height = unit(4, "cm")
  )
)

draw(ht, heatmap_legend_side = "right")


