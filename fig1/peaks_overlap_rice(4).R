library(GenomicRanges)
library(rtracklayer)
library(ChIPpeakAnno)
library(ggVennDiagram)
library(scales)  # for alpha()
library(ggplot2)

peak_1_rep1 <- read.table(
  "D:/sn-CT_broadpeaks/5/NA_peaks.broadPeak",
  col.names = c("chrom", "start", "end", "name", "score", "strand", 
                "signal", "pval", "qval")
)

peak_1_rep2 <- read.table(
  "D:/sn-CT_broadpeaks/6/NA_peaks.broadPeak",
  col.names = c("chrom", "start", "end", "name", "score", "strand", 
                "signal", "pval", "qval")
)

peak_2_rep1 <- read.table(
  "D:/sn-CT_broadpeaks/7/NA_peaks.broadPeak",
  col.names = c("chrom", "start", "end", "name", "score", "strand", 
                "signal", "pval", "qval")
)

peak_2_rep2 <- read.table(
  "D:/sn-CT_broadpeaks/8/NA_peaks.broadPeak",
  col.names = c("chrom", "start", "end", "name", "score", "strand", 
                "signal", "pval", "qval")
)

peak_3_rep1 <- read.table(
  "C:/Users/admin/Desktop/chip-peak/1611.broadPeak",
  col.names = c("chrom", "start", "end", "name", "score", "strand", 
                "signal", "pval", "qval")
)

peak_3_rep2 <- read.table(
  "C:/Users/admin/Desktop/chip-peak/1612.broadPeak",
  col.names = c("chrom", "start", "end", "name", "score", "strand", 
                "signal", "pval", "qval")
)

peak_4_rep1 <- read.table(
  "C:/Users/admin/Desktop/chip-peak/1645.broadPeak",
  col.names = c("chrom", "start", "end", "name", "score", "strand", 
                "signal", "pval", "qval")
)

peak_4_rep2 <- read.table(
  "C:/Users/admin/Desktop/chip-peak/1646.broadPeak",
  col.names = c("chrom", "start", "end", "name", "score", "strand", 
                "signal", "pval", "qval")
)


gr.peaks_1_1 <- makeGRangesFromDataFrame(peak_1_rep1)
gr.peaks_1_2 <- makeGRangesFromDataFrame(peak_1_rep2)
gr.peaks_2_1 <- makeGRangesFromDataFrame(peak_2_rep1)
gr.peaks_2_2 <- makeGRangesFromDataFrame(peak_2_rep2)
gr.peaks_3_1 <- makeGRangesFromDataFrame(peak_3_rep1)
gr.peaks_3_2 <- makeGRangesFromDataFrame(peak_3_rep2)
gr.peaks_4_1 <- makeGRangesFromDataFrame(peak_4_rep1)
gr.peaks_4_2 <- makeGRangesFromDataFrame(peak_4_rep2)

peaks_group1 <- reduce(x = c(gr.peaks_1_1, gr.peaks_1_2))
peaks_group2 <- reduce(x = c(gr.peaks_2_1, gr.peaks_2_2))
peaks_group3 <- reduce(x = c(gr.peaks_3_1, gr.peaks_3_2))
peaks_group4 <- reduce(x = c(gr.peaks_4_1, gr.peaks_4_2))


overlap_result <- findOverlapsOfPeaks(
  group1 = peaks_group1,
  group2 = peaks_group2,
  group3 = peaks_group3,
  group4 = peaks_group4
)


pdf("venn_diagram_peaks_rice(4).pdf", width = 15, height = 15)
colors <- c("#264653","#2A9D8F","#E9C46A","#F4A261")
sample_names <- c("snCT_Nano-Tn5", "snCT_Nano-Tn5_Beads", "Bulk_ChIP_yl", "Bulk_ChIP_ml")
makeVennDiagram(
  overlap_result,
  NameOfPeaks = sample_names,
  fill = colors,
  lty = "blank",
  main = "Overlap of Peaks",
  cat.col = colors,           # 设置分类标签的颜色
  cat.cex = 1.2,              # 分类标签字体大小
  cex = 1.2,                  # 数字标签字体大小
  main.cex = 1.5              # 标题字体大小
) +
  scale_fill_gradientn(colours = alpha(colors,0.3))

dev.off()












