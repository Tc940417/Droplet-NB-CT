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
library(Seurat)
nip <- readRDS("D:/sn-CUT&TAG-BGI_result/1_scRNA/1_scRNA/nip.gene_500-5000.resolution_0.25_anno_unknown.rds")
nip$label<- paste(nip$cell_type, nip$group, sep = "-")
nip <- JoinLayers(nip)
Idents(nip) <- nip@meta.data$label

Me_degs <- FindMarkers(nip, 
                       ident.1 = "Mesophyll-case", 
                       ident.2 = "Mesophyll-ctrl", 
                       test.use="MAST", 
                       min.pct= 0.1, 
                       logfc.threshold= 0.25,
                       latent.vars = "nCount_RNA",
                       slot = "data")
Lp_degs <- FindMarkers(nip, 
                       ident.1 = "Large parenchyma-case", 
                       ident.2 = "Large parenchyma-ctrl", 
                       test.use="MAST", 
                       min.pct= 0.1, 
                       logfc.threshold= 0.25,
                       latent.vars = "nCount_RNA",
                       slot = "data")
Ep_degs <- FindMarkers(nip, 
                       ident.1 = "Epidermis-case", 
                       ident.2 = "Epidermis-ctrl", 
                       test.use="MAST", 
                       min.pct= 0.1, 
                       logfc.threshold= 0.25,
                       latent.vars = "nCount_RNA",
                       slot = "data")
Vs_degs <- FindMarkers(nip, 
                       ident.1 = "Vascular-case", 
                       ident.2 = "Vascular-ctrl", 
                       test.use="MAST", 
                       min.pct= 0.1, 
                       logfc.threshold= 0.25,
                       latent.vars = "nCount_RNA",
                       slot = "data")
Gc_degs <- FindMarkers(nip, 
                       ident.1 = "Guard cell-case", 
                       ident.2 = "Guard cell-ctrl", 
                       test.use="MAST", 
                       min.pct= 0.1, 
                       logfc.threshold= 0.25,
                       latent.vars = "nCount_RNA",
                       slot = "data")
Fb_degs <- FindMarkers(nip, 
                       ident.1 = "Fiber-case", 
                       ident.2 = "Fiber-ctrl", 
                       test.use="MAST", 
                       min.pct= 0.1, 
                       logfc.threshold= 0.25,
                       latent.vars = "nCount_RNA",
                       slot = "data")

write.csv(Me_degs,
          file = "Mesophyll_case_vs_ctrl.csv")
write.csv(Lp_degs,
          file = "LargeParenchyma_case_vs_ctrl.csv")
write.csv(Ep_degs,
          file = "Epidermis_case_vs_ctrl.csv")
write.csv(Vs_degs,
          file = "Vascular_case_vs_ctrl.csv")
write.csv(Gc_degs,
          file = "GuardCell_case_vs_ctrl.csv")
write.csv(Fb_degs,
          file = "Fiber_case_vs_ctrl.csv")

Me_degs.sig <- Me_degs[Me_degs$p_val_adj < 0.05, ]
Lp_degs.sig <- Lp_degs[Lp_degs$p_val_adj < 0.05, ]
Ep_degs.sig <- Ep_degs[Ep_degs$p_val_adj < 0.05, ]
Vs_degs.sig <- Vs_degs[Vs_degs$p_val_adj < 0.05, ]
Gc_degs.sig <- Gc_degs[Gc_degs$p_val_adj < 0.05, ]
Fb_degs.sig <- Fb_degs[Fb_degs$p_val_adj < 0.05, ]

deg.list <- list(
  Mesophyll       = Me_degs,
  LargeParenchyma = Lp_degs,
  Epidermis       = Ep_degs,
  Vascular        = Vs_degs,
  GuardCell       = Gc_degs,
  Fiber           = Fb_degs
)
log2FC <- 1

deg.split <- lapply(deg.list, function(df) {
  list(
    up = df[df$p_val_adj < 0.05 & df$avg_log2FC >=  log2FC, ],
    down = df[df$p_val_adj < 0.05 & df$avg_log2FC <= -log2FC, ],
    ns = df[!(df$p_val_adj < 0.05 & abs(df$avg_log2FC) >= log2FC), ]
  )
})
sapply(deg.split, function(x) {
  c(
    up   = nrow(x$up),
    down = nrow(x$down),
    ns   = nrow(x$ns)
  )
})

out.dir <- "E:/sn_CT/table/rice_snRNA_DEG/DEGs_split"

for (ct in names(deg.split)) {
  ct.dir <- file.path(out.dir, ct)
  dir.create(ct.dir, showWarnings = FALSE)
  
  write.csv(deg.split[[ct]]$up,
            file = file.path(ct.dir, "up.csv"),
            row.names = TRUE)
  
  write.csv(deg.split[[ct]]$down,
            file = file.path(ct.dir, "down.csv"),
            row.names = TRUE)
  
  write.csv(deg.split[[ct]]$ns,
            file = file.path(ct.dir, "ns.csv"),
            row.names = TRUE)
}























