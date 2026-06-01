library (Seurat)
library (harmony)
library (ggplot2)
library (viridis)
library (cowplot)
library (dplyr)
library (tidydr)
library (stringr)
library (scCustomize)


cnt <- Read10X(data.dir="../3_clean_mtx/NIP_TR197_none/", gene.column=1)
ctrl1 <- CreateSeuratObject (counts=cnt, min.cells=1, min.features=1)

cnt <- Read10X(data.dir="../3_clean_mtx/NIP_TR198_none/", gene.column=1)
ctrl2 <- CreateSeuratObject (counts=cnt, min.cells=1, min.features=1)

cnt <- Read10X(data.dir="../3_clean_mtx/NIP_TR199_45C90min/", gene.column=1)
case1 <- CreateSeuratObject (counts=cnt, min.cells=1, min.features=1)

cnt <- Read10X(data.dir="../3_clean_mtx/NIP_TR200_45C90min/", gene.column=1)
case2 <- CreateSeuratObject (counts=cnt, min.cells=1, min.features=1)

ctrl1@meta.data$spl <- "ctrl1"
ctrl2@meta.data$spl <- "ctrl2"
case1@meta.data$spl <- "case1"
case2@meta.data$spl <- "case2"

ctrl1@meta.data$group <- "ctrl"
ctrl2@meta.data$group <- "ctrl"
case1@meta.data$group <- "case"
case2@meta.data$group <- "case"

pdf ("nip.gene_500-5000.resolution_0.5.pdf");
obj <- merge (ctrl1, y=c(ctrl2, case1, case2), add.cell.ids=c("ctrl1", "ctrl2", "case1", "case2"))
obj <- subset (obj, subset= nFeature_RNA>500 & nFeature_RNA<8000);

obj <- NormalizeData (obj);
obj <- FindVariableFeatures (obj);
obj <- ScaleData (obj);
obj <- RunPCA (obj);

obj <- IntegrateLayers (obj, method=CCAIntegration,
		                    orig.reduction="pca", new.reduction = "integrated.cca");
obj <- FindNeighbors (obj, reduction="integrated.cca", dims=1:30);
obj <- FindClusters (obj, resolution=0.25);
obj <- RunUMAP (obj, reduction="integrated.cca", dims=1:30);

saveRDS (obj, file="nip.gene_500-5000.resolution_0.5.rds");
saveRDS (nip.anno, file="nip.gene_500-5000.resolution_0.5.anno.rds");
nip <-readRDS("D:/sn-CUT&TAG-BGI_result/1_scRNA/1_scRNA/nip.gene_500-5000.resolution_0.25_anno.rds")
nip <-readRDS("D:/sn-CUT&TAG-BGI_result/1_scRNA/1_scRNA/nip.gene_500-5000.resolution_0.25_anno_unknown.rds")
DimPlot (obj, group.by="group") + ggtitle(label="CCA Integration UMAP");
DimPlot (obj, group.by="spl") + ggtitle(label="CCA Integration UMAP");
DimPlot (nip, label=T) + ggtitle(label="CCA Integration UMAP");

nip <- JoinLayers (nip);
deg <- FindAllMarkers (nip, only.pos=F);
write.csv(
  deg,
  file = "E:/sn_CT/table/nip..0.25.csv",
  row.names = TRUE,
  quote = FALSE
)
dev.off ()

