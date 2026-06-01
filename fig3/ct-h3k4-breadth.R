library(ggplot2)
library(tidyr)
library(Seurat)
library(Signac)
library(patchwork)
library(dplyr)
library(purrr)
library(dplyr)
library(igraph)
library(data.table)
RNA <- readRDS("D:/sn-CUT&TAG-BGI_result/1_scRNA/1_scRNA/nip.gene_500-5000.resolution_0.25_anno_unknown.rds")
ct_RNA <- subset(RNA , subset = group == "ctrl")
ct_RNA <- JoinLayers(ct_RNA)
deg <- FindAllMarkers(ct_RNA, test.use="MAST", min.pct= 0.1, logfc.threshold= 0.5,latent.vars = "nCount_RNA", only.pos = TRUE)
write.csv(deg,
          file = "E:/sn_CT/table/rice_snRNA_DEG/deg_FindAllMarkers_MAST.csv",
          row.names = FALSE)
deg_filt <- deg %>%
  filter(avg_log2FC > 0.5, p_val_adj < 0.05) %>%
  mutate(gene = gsub("-", "_", gene))



h3k4 <- readRDS("E:/sn_CT/RDS&h5ad/snCTaddRNAaddChromvar/ct.h3k4.rds")
DefaultAssay(h3k4) <- "ATAC"
library(Signac)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
gene_anno <- Annotation(h3k4)
gene_anno <- gene_anno[gene_anno$type == "gene"]
promoters <- promoters(
  gene_anno,
  upstream = 2000,
  downstream = 2000
)
peaks <- granges(h3k4)
peaks$breadth <- width(peaks)
overlap <- findOverlaps(promoters, peaks)

peak_gene_df <- data.frame(
  gene = promoters$gene_name[queryHits(overlap)],
  breadth = peaks$breadth[subjectHits(overlap)]
)

gene_breadth <- peak_gene_df %>%
  group_by(gene) %>%
  summarise(
    peak_breadth = sum(breadth)
  )

all_genes_breadth <- gene_breadth %>%
  mutate(group = "All genes")

deg_genes_breadth <- gene_breadth %>%
  filter(gene %in% deg_filt$gene) %>%
  mutate(group = "DEG genes")

plot_df <- bind_rows(all_genes_breadth, deg_genes_breadth)


ggplot(plot_df, aes(x = peak_breadth, fill = group)) +
  geom_histogram(
    aes(y = ..density..),
    bins = 50,
    position = "identity",
    alpha = 0.55,
    color = "#3f3f3f",   # 柱子线框
    linewidth = 0.3
  ) +
  scale_x_log10() +
  theme_classic(base_size = 12) +
  labs(
    x = "H3K4me3 peak breadth (bp)",
    y = "Density",
    fill = ""
  ) +
  scale_fill_manual(
    values = c(
      "All genes" = "#c97c6d",   # 莫兰迪豆沙红
      "DEG genes" = "#7f9d8c"    # 莫兰迪灰绿
    )
  )+ geom_density(
    aes(color = group),
    size = 0.7,
    fill = NA
  ) +
  scale_color_manual(
    values = c(
      "All genes" = "#9b5f54",
      "DEG genes" = "#5f7f6b"
    )
  )






deg_anno <- deg_filt %>%
  dplyr::select(gene, cluster) %>%
  distinct()

deg_genes_breadth <- gene_breadth %>%
  inner_join(deg_anno, by = "gene") %>%
  mutate(group = cluster)


all_genes_breadth <- gene_breadth %>%
  mutate(group = "All genes")


plot_df <- bind_rows(all_genes_breadth, deg_genes_breadth)


distinguishable_morandi_colors <- c(
  "Guard cell"        = "#d69971",
  "Mesophyll"         = "#64a776",
  "Epidermis"         = "#b85292",
  "Large parenchyma"  = "#537eb7",
  "Fiber"             = "#efd2c9",
  "Vascular"          = "#9f8d89",
  "All genes"         = "#000000"  # 黑色
)

plot_df$group <- factor(
  plot_df$group,
  levels = c(
    "Large parenchyma",
    "Mesophyll",
    "Epidermis",
    "Vascular",
    "Fiber",
    "Guard cell",
    "All genes"
  )
)

ggplot() +
  stat_ecdf(
    data = subset(plot_df, group != "All genes"),
    aes(x = peak_breadth, color = group),
    size = 0.8
  ) +
  stat_ecdf(
    data = subset(plot_df, group == "All genes"),
    aes(x = peak_breadth, color = group),
    size = 0.5,
    linetype = "dashed"
  ) +
  scale_x_continuous(
    limits = c(0, 5000),
    breaks = c(0, 1000, 2000, 3000, 4000, 5000)
  ) +
  scale_color_manual(values = distinguishable_morandi_colors) +
  theme_classic(base_size = 12) +
  labs(
    x = "H3K4me3 peak breadth (bp)",
    y = "Cumulative fraction",
    color = "Cell type"
  )



