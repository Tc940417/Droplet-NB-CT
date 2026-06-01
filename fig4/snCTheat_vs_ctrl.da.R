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
h3k4 <- readRDS("E:/sn_CT/RDS&h5ad/snCTaddRNAaddChromvar/h3k4.rds")
h3k27 <- readRDS("E:/sn_CT/RDS&h5ad/snCTaddRNAaddChromvar/h3k27.rds")
DefaultAssay(h3k4) <- "ATAC"
DefaultAssay(h3k27) <- "ATAC"

h3k4 <- RunTFIDF(h3k4)
h3k27 <- RunTFIDF(h3k27)

h3k4$celltype <- paste(h3k4$batch, h3k4$cell_type, sep = "_")
h3k27$celltype <- paste(h3k27$batch, h3k27$cell_type, sep = "_")

Idents(h3k4) <- "celltype"
Idents(h3k27) <- "celltype"



k4.me <- FindMarkers(
  object = h3k4,
  ident.1 = "heat_h3k4me3_Mesophyll",
  ident.2 = "ctrl_h3k4me3_Mesophyll",
  test.use = 'LR',
  min.pct = 0.1,
  latent.vars = c("nCount_ATAC","nFeature_ATAC"),
  logfc.threshold= 0.25,
  slot = "data"
)

k4.lp <- FindMarkers(
  object = h3k4,
  ident.1 = "heat_h3k4me3_Large parenchyma",
  ident.2 = "ctrl_h3k4me3_Large parenchyma",
  test.use = 'LR',
  min.pct = 0.1,
  latent.vars = c("nCount_ATAC","nFeature_ATAC"),
  logfc.threshold= 0.25,
  slot = "data"
)

k4.gc <- FindMarkers(
  object = h3k4,
  ident.1 = "heat_h3k4me3_Guard cell",
  ident.2 = "ctrl_h3k4me3_Guard cell",
  test.use = 'LR',
  min.pct = 0.1,
  latent.vars = c("nCount_ATAC","nFeature_ATAC"),
  logfc.threshold= 0.25,
  slot = "data"
)

k4.vs <- FindMarkers(
  object = h3k4,
  ident.1 = "heat_h3k4me3_Vascular",
  ident.2 = "ctrl_h3k4me3_Vascular",
  test.use = 'LR',
  min.pct = 0.1,
  latent.vars = c("nCount_ATAC","nFeature_ATAC"),
  logfc.threshold= 0.25,
  slot = "data"
)

k4.ep <- FindMarkers(
  object = h3k4,
  ident.1 = "heat_h3k4me3_Epidermis",
  ident.2 = "ctrl_h3k4me3_Epidermis",
  test.use = 'LR',
  min.pct = 0.1,
  latent.vars = c("nCount_ATAC","nFeature_ATAC"),
  logfc.threshold= 0.25,
  slot = "data"
)

k4.fb <- FindMarkers(
  object = h3k4,
  ident.1 = "heat_h3k4me3_Fiber",
  ident.2 = "ctrl_h3k4me3_Fiber",
  test.use = 'LR',
  min.pct = 0.1,
  latent.vars = c("nCount_ATAC","nFeature_ATAC"),
  logfc.threshold= 0.25,
  slot = "data"
)

########################################################

k27.me <- FindMarkers(
  object = h3k27,
  ident.1 = "heat_h3k27ac_Mesophyll",
  ident.2 = "ctrl_h3k27ac_Mesophyll",
  test.use = 'LR',
  min.pct = 0.1,
  latent.vars = c("nCount_ATAC","nFeature_ATAC"),
  logfc.threshold= 0.25,
  slot = "data"
)

k27.lp <- FindMarkers(
  object = h3k27,
  ident.1 = "heat_h3k27ac_Large parenchyma",
  ident.2 = "ctrl_h3k27ac_Large parenchyma",
  test.use = 'LR',
  min.pct = 0.1,
  latent.vars = c("nCount_ATAC","nFeature_ATAC"),
  logfc.threshold= 0.25,
  slot = "data"
)

k27.gc <- FindMarkers(
  object = h3k27,
  ident.1 = "heat_h3k27ac_Guard cell",
  ident.2 = "ctrl_h3k27ac_Guard cell",
  test.use = 'LR',
  min.pct = 0.1,
  latent.vars = c("nCount_ATAC","nFeature_ATAC"),
  logfc.threshold= 0.25,
  slot = "data"
)

k27.vs <- FindMarkers(
  object = h3k27,
  ident.1 = "heat_h3k27ac_Vascular",
  ident.2 = "ctrl_h3k27ac_Vascular",
  test.use = 'LR',
  min.pct = 0.1,
  latent.vars = c("nCount_ATAC","nFeature_ATAC"),
  logfc.threshold= 0.25,
  slot = "data"
)

k27.ep <- FindMarkers(
  object = h3k27,
  ident.1 = "heat_h3k27ac_Epidermis",
  ident.2 = "ctrl_h3k27ac_Epidermis",
  test.use = 'LR',
  min.pct = 0.1,
  latent.vars = c("nCount_ATAC","nFeature_ATAC"),
  logfc.threshold= 0.25,
  slot = "data"
)

k27.fb <- FindMarkers(
  object = h3k27,
  ident.1 = "heat_h3k27ac_Fiber",
  ident.2 = "ctrl_h3k27ac_Fiber",
  test.use = 'LR',
  min.pct = 0.1,
  latent.vars = c("nCount_ATAC","nFeature_ATAC"),
  logfc.threshold= 0.25,
  slot = "data"
)



library(dplyr)
library(tibble)

k27.list <- list(
  Mesophyll = k27.me,
  Large_parenchyma = k27.lp,
  Guard_cell = k27.gc,
  Vascular = k27.vs,
  Epidermis = k27.ep,
  Fiber = k27.fb
)

k27.sig <- lapply(k27.list, function(df) {
  df %>%
    rownames_to_column("peak") %>%
    filter(
      abs(avg_log2FC) > 1,
      p_val_adj < 0.05,
      pct.1 >= 0.1
    )
})

k27.sig.df <- bind_rows(
  lapply(names(k27.sig), function(ct) {
    k27.sig[[ct]] %>%
      mutate(celltype = ct)
  })
)

write.csv(k4.me, "marker_results/H3K4me3/k4_Mesophyll.csv")
write.csv(k4.lp, "marker_results/H3K4me3/k4_LargeParenchyma.csv")
write.csv(k4.gc, "marker_results/H3K4me3/k4_GuardCell.csv")
write.csv(k4.vs, "marker_results/H3K4me3/k4_Vascular.csv")
write.csv(k4.ep, "marker_results/H3K4me3/k4_Epidermis.csv")
write.csv(k4.fb, "marker_results/H3K4me3/k4_Fiber.csv")

write.csv(k27.me, "marker_results/H3K27ac/k27_Mesophyll.csv")
write.csv(k27.lp, "marker_results/H3K27ac/k27_LargeParenchyma.csv")
write.csv(k27.gc, "marker_results/H3K27ac/k27_GuardCell.csv")
write.csv(k27.vs, "marker_results/H3K27ac/k27_Vascular.csv")
write.csv(k27.ep, "marker_results/H3K27ac/k27_Epidermis.csv")
write.csv(k27.fb, "marker_results/H3K27ac/k27_Fiber.csv")




