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
h3k4 <- readRDS("E:/sn_CT/RDS&h5ad/snCTaddRNAaddChromvar/h3k4.rds")
DefaultAssay(h3k4) <- "ATAC"
h3k4 <- RunTFIDF(h3k4)
h3k4$batch[h3k4$batch == "ctrl_h3k4me3"] <- "ctrl"
h3k4$batch[h3k4$batch == "heat_h3k4me3"] <- "heat"
h3k4$label<- paste(h3k4$cell_type, h3k4$batch, sep = "-")
Idents(h3k4) <- h3k4@meta.data$label

pb <- AggregateExpression(
  object = h3k4,
  assays = "ATAC",
  group.by = c("cell_type", "batch", "dataset"),
  slot = "counts",
  return.seurat = TRUE
)

count.mat <- GetAssayData(pb, assay = "ATAC", slot = "counts")
coldata <- pb@meta.data[, c("cell_type", "batch", "dataset")]

##############################################################################
ct <- "Vascular"
keep <- coldata$cell_type == ct
count.ct <- count.mat[, keep]
coldata.ct <- coldata[keep, ]
coldata.ct$batch <- factor(coldata.ct$batch, levels = c("ctrl", "heat"))
coldata.ct$dataset <- factor(coldata.ct$dataset)

library(DESeq2)
dds <- DESeqDataSetFromMatrix(
  countData = count.ct,
  colData   = coldata.ct,
  design    = ~ batch
)
dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)
res <- as.data.frame(
  results(dds, contrast = c("batch", "heat", "ctrl"))
)


##############################################################################



da <- FindMarkers(
  h3k4,
  ident.1 = "Vascular-heat",
  ident.2 = "Vascular-ctrl",
  test.use = "LR",
  latent.vars = c("nCount_ATAC", "nFeature_ATAC"),
  min.pct = 0.05,
  logfc.threshold = 0,
  only.pos = FALSE,
  slot = "data"
)










