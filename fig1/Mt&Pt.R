library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(clustree)
library(rtracklayer)

# 文件夹编号
folders <- c(5:8)

# 列名
colnames_broadpeak <- c("chrom", "start", "end", "name", "score", "strand", 
                        "signal", "pval", "qval")

for (i in folders) {
  # 构建文件路径
  file_path <- file.path("D:/sn-CT_broadpeaks", i, "NA_peaks.broadPeak")
  
  # 读取表格
  peaks <- read.table(file_path, col.names = colnames_broadpeak)
  
  # 构建 GRanges
  gr <- makeGRangesFromDataFrame(peaks, keep.extra.columns = TRUE)
  
  # 动态赋值为 peaks_3, peaks_4, ...
  assign(paste0("peaks_", i), peaks)
  assign(paste0("gr.peaks_", i), gr)
}

for (i in c(5:8)) {
  gr_obj <- get(paste0("gr.peaks_", i))                # 取出对象
  peakwidths <- width(gr_obj)                          # 计算宽度
  gr_obj <- gr_obj[peakwidths < 10000 & peakwidths > 20]  # 筛选长度
  assign(paste0("gr.peaks_", i), gr_obj)               # 赋回原变量
}



for (i in folders) {
  file_path <- file.path("D:/sn-CT_broadpeaks", i, "singlecell.csv")
  
  # 读取 CSV，设置行名
  md <- read.table(file = file_path,
                   stringsAsFactors = FALSE,
                   sep = ",",
                   header = TRUE,
                   row.names = 1)
  
  # 筛选 is_cell_barcode > 0 的行
  md <- md[md$is_cell_barcode > 0, ]
  
  # 赋值为 md_3、md_4 等变量名
  assign(paste0("md_", i), md)
}

for (i in c(5:8)){
  md <- get(paste0("md_", i))  # 获取对象
  
  # 添加三列
  md$log10_uniqueFrags <- log10(md$fragments)
  md$pct_reads_in_peaks <- md$peak_region_fragments / md$fragments * 100
  md$pct_reads_in_tss <- md$TSS_region_fragments / md$fragments * 100
  
  # 更新变量
  assign(paste0("md_", i), md)
}




for (i in c(5:8)) {
  # 构建文件路径
  frag_path <- file.path("D:/sn-CT_broadpeaks", i, "include.MtPt.fragments.tsv.gz")
  
  # 获取细胞条码
  md <- get(paste0("md_", i))
  barcodes <- rownames(md)
  
  # 创建 FragmentObject
  frag_obj <- CreateFragmentObject(
    path = frag_path,
    cells = barcodes
  )
  
  # 赋值为 frags_3、frags_4 等
  assign(paste0("frags_", i), frag_obj)
}


for (i in c(5:8)) {
  message("Processing sample ", i)
  
  # 获取对象
  frag <- get(paste0("frags_", i))
  peaks <- get(paste0("gr.peaks_", i))
  cells <- rownames(get(paste0("md_", i)))
  
  # 创建 count matrix
  counts <- FeatureMatrix(
    fragments = frag,
    features = peaks,
    cells = cells
  )
  
  # 保存为 counts_i
  assign(paste0("counts_", i), counts)
}

for (i in c(5:8)) {
  message("Creating Seurat object for sample ", i)
  
  # 提取变量名
  counts <- get(paste0("counts_", i))
  frag <- get(paste0("frags_", i))
  md <- get(paste0("md_", i))
  
  # 创建 ChromatinAssay
  assay <- CreateChromatinAssay(
    counts = counts,
    fragments = frag
  )
  
  # 创建 Seurat 对象
  seurat_obj <- CreateSeuratObject(
    counts = assay,
    assay = "CT",
    meta.data = md
  )
  
  # 保存对象
  assign(paste0("assay_", i), assay)
  assign(paste0("snCT_", i), seurat_obj)
}

snCT_5[["percent.mt"]] <- PercentageFeatureSet(snCT_5, pattern = "^Mt")
snCT_5[["percent.pt"]] <- PercentageFeatureSet(snCT_5, pattern = "^Pt")

snCT_6[["percent.mt"]] <- PercentageFeatureSet(snCT_6, pattern = "^Mt")
snCT_6[["percent.pt"]] <- PercentageFeatureSet(snCT_6, pattern = "^Pt")

snCT_7[["percent.mt"]] <- PercentageFeatureSet(snCT_7, pattern = "^Mt")
snCT_7[["percent.pt"]] <- PercentageFeatureSet(snCT_7, pattern = "^Pt")

snCT_8[["percent.mt"]] <- PercentageFeatureSet(snCT_8, pattern = "^Mt")
snCT_8[["percent.pt"]] <- PercentageFeatureSet(snCT_8, pattern = "^Pt")

snCT_13[["percent.mt"]] <- PercentageFeatureSet(snCT_13, pattern = "^Mt")
snCT_13[["percent.pt"]] <- PercentageFeatureSet(snCT_13, pattern = "^Pt")

snCT_14[["percent.mt"]] <- PercentageFeatureSet(snCT_14, pattern = "^Mt")
snCT_14[["percent.pt"]] <- PercentageFeatureSet(snCT_14, pattern = "^Pt")

snCT_15[["percent.mt"]] <- PercentageFeatureSet(snCT_15, pattern = "^Mt")
snCT_15[["percent.pt"]] <- PercentageFeatureSet(snCT_15, pattern = "^Pt")

snCT_16[["percent.mt"]] <- PercentageFeatureSet(snCT_16, pattern = "^Mt")
snCT_16[["percent.pt"]] <- PercentageFeatureSet(snCT_16, pattern = "^Pt")

snCT_17[["percent.mt"]] <- PercentageFeatureSet(snCT_17, pattern = "^Mt")
snCT_17[["percent.pt"]] <- PercentageFeatureSet(snCT_17, pattern = "^Pt")

snCT_18[["percent.mt"]] <- PercentageFeatureSet(snCT_18, pattern = "^Mt")
snCT_18[["percent.pt"]] <- PercentageFeatureSet(snCT_18, pattern = "^Pt")


##使用GTF文件创建注释文件##
annotations <- import("D:/sn-CT/Osref/Oryza_sativa.MSU.chr.checked.signac.gtf")

##设置染色体长度##
seqlengths <- c("1"=43270923, "2"=35937250, "3"=36413819, "4"=35502694, "5"=29958434, "6"=31248787, "7"=29697621,
                "8"=28443022 ,"9"=23012720, "10"=23207287, "11"=29021106, "12"=27531856) 
seqlengths(annotations) <- seqlengths

##增加gene_biotype列，默认全部设为蛋白编码##
annotations$gene_biotype <- "protein_coding"
genome(annotations) <- "MSU_osa1r7" 
Annotation(snCT_5) <- annotations
Annotation(snCT_6) <- annotations
Annotation(snCT_7) <- annotations
Annotation(snCT_8) <- annotations


for (i in c(5:8)) {
  message("Processing snCT_", i)
  
  obj <- get(paste0("snCT_", i))  # 正确变量名
  
  obj <- NucleosomeSignal(object = obj)
  obj <- TSSEnrichment(object = obj)
  
  obj$nucleosome_group <- ifelse(
    obj$nucleosome_signal > 4, 'NS > 4', 'NS < 4'
  )
  
  assign(paste0("snCT_", i), obj)  # 正确赋值
}


for (i in c(5:8)) {
  message("Filtering snCT_", i)
  
  varname <- paste0("snCT_", i)
  
  if (!exists(varname)) {
    warning("Object ", varname, " not found, skipping.")
    next
  }
  
  obj <- get(varname)
  
  obj <- subset(
    x = obj,
    subset = nCount_CT > 1000 &
      nCount_CT < 100000 &
      pct_reads_in_peaks > 15 &
      nucleosome_signal < 4 &
      TSS.enrichment > 1
  )
  
  assign(varname, obj)
}

snCT_5 <- FRiP(object = snCT_5, assay = "CT", total.fragments = "fragments", col.name = "FRiP", verbose = TRUE)
snCT_6 <- FRiP(object = snCT_6, assay = "CT", total.fragments = "fragments", col.name = "FRiP", verbose = TRUE)
snCT_7 <- FRiP(object = snCT_7, assay = "CT", total.fragments = "fragments", col.name = "FRiP", verbose = TRUE)
snCT_8 <- FRiP(object = snCT_8, assay = "CT", total.fragments = "fragments", col.name = "FRiP", verbose = TRUE)

saveRDS(snCT_5, file = "E:/sn_CT/RDS&h5ad/1.qc/nano-1.rds")
saveRDS(snCT_6, file = "E:/sn_CT/RDS&h5ad/1.qc/nano-2.rds")
saveRDS(snCT_7, file = "E:/sn_CT/RDS&h5ad/1.qc/NB-1.rds")
saveRDS(snCT_8, file = "E:/sn_CT/RDS&h5ad/1.qc/NB-2.rds")



seurat_list <- list(sample1_rep1 = snCT_5,sample1_rep2 = snCT_6,
                    sample2_rep1 = snCT_7,sample2_rep2 = snCT_8)
##Mt##

mt_df <- imap_dfr(seurat_list, ~{
  data.frame(
    cell_id = colnames(.x),
    mt = .x$percent.mt,
    sample = .y,
    stringsAsFactors = FALSE
  )
}) %>% 
  tidyr::pivot_wider(
    names_from = sample,
    values_from = mt
  )


mt_df_1 <- mt_df[,-1]

colnames(mt_df_1) <- c("ctrl_H3K4me3_Nano_Tn5_rep1", "ctrl_H3K4me3_Nano_Tn5_rep2", 
                       "ctrl_H3K4me3_Nano_Tn5_Beads_rep1","ctrl_H3K4me3_Nano_Tn5_Beads_rep2", 
                      "ctrl_H3K27ac_Nano_Tn5_Beads_rep1", "ctrl_H3K27ac_Nano_Tn5_Beads_rep2", 
                      "heat_H3K4me3_Nano_Tn5_Beads_rep1", "heat_H3K4me3_Nano_Tn5_Beads_rep2", 
                      "heat_H3K27ac_Nano_Tn5_Beads_rep1", "heat_H3K27ac_Nano_Tn5_Beads_rep2" 
)

sample_colors <- c(
  rep('#df5734', 2),  # 
  rep('#64a776', 2),  # 
  rep('#2d3462', 2),  # 
  rep('#ebb1a4', 2),  # 
  rep('#76a2be', 2)   # 
)

plot_long <- mt_df_1 %>%
  pivot_longer(
    cols = everything(),
    names_to = "Sample",
    values_to = "percent.mt"
  ) %>%
  mutate(
    Sample = factor(Sample, levels = colnames(mt_df_1)), # 保持原始列顺序
    Group = rep(1:5, each = 2)[as.numeric(Sample)] # 创建分组变量用于颜色映射
  )

plot <- ggplot(plot_long, aes(x = Sample, y = percent.mt , fill = factor(Group))) +
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
  # 美化主题
  labs(
    title = "Mitochondrial percentage",
    x = NULL,
    y = NULL,
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
    yintercept = c(1, 3, 5), 
    linetype = "dashed", 
    color = "gray50",
    alpha = 0.5
  )
plot
ggsave(
  filename = "rice_Mt.percent_violin_plot.pdf",
  plot = plot,
  device = "pdf",
  width = 10,        # 宽度（英寸）
  height = 6,        # 高度（英寸）
  units = "in",      # 尺寸单位
  dpi = 1200          # 分辨率
)
















##Pt##

pt_df <- imap_dfr(seurat_list, ~{
  data.frame(
    cell_id = colnames(.x),
    pt = .x$percent.pt,
    sample = .y,
    stringsAsFactors = FALSE
  )
}) %>% 
  tidyr::pivot_wider(
    names_from = sample,
    values_from = pt
  )


pt_df_1 <- pt_df[,-1]

colnames(pt_df_1) <- c("ctrl_H3K4me3_Nano_Tn5_rep1", "ctrl_H3K4me3_Nano_Tn5_rep2", 
                       "ctrl_H3K4me3_Nano_Tn5_Beads_rep1","ctrl_H3K4me3_Nano_Tn5_Beads_rep2", 
                       "ctrl_H3K27ac_Nano_Tn5_Beads_rep1", "ctrl_H3K27ac_Nano_Tn5_Beads_rep2", 
                       "heat_H3K4me3_Nano_Tn5_Beads_rep1", "heat_H3K4me3_Nano_Tn5_Beads_rep2", 
                       "heat_H3K27ac_Nano_Tn5_Beads_rep1", "heat_H3K27ac_Nano_Tn5_Beads_rep2" 
)

sample_colors <- c(
  rep('#df5734', 2),  # 
  rep('#64a776', 2),  # 
  rep('#2d3462', 2),  # 
  rep('#ebb1a4', 2),  # 
  rep('#76a2be', 2)   # 
)

plot_long <- pt_df_1 %>%
  pivot_longer(
    cols = everything(),
    names_to = "Sample",
    values_to = "percent.pt"
  ) %>%
  mutate(
    Sample = factor(Sample, levels = colnames(pt_df_1)), # 保持原始列顺序
    Group = rep(1:5, each = 2)[as.numeric(Sample)] # 创建分组变量用于颜色映射
  )

plot <- ggplot(plot_long, aes(x = Sample, y = percent.pt , fill = factor(Group))) +
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
  # 美化主题
  labs(
    title = "Chloroplast percentage",
    x = NULL,
    y = NULL,
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
    yintercept = c(0.5, 1, 1.5), 
    linetype = "dashed", 
    color = "gray50",
    alpha = 0.5
  )
plot
ggsave(
  filename = "rice_Pt.percent_violin_plot.pdf",
  plot = plot,
  device = "pdf",
  width = 10,        # 宽度（英寸）
  height = 6,        # 高度（英寸）
  units = "in",      # 尺寸单位
  dpi = 1200          # 分辨率
)


































