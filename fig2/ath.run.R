library(Seurat)

cnt91 <- Read10X (data.dir="../3_clean_mtx/Ara_TR091_none/", gene.column=1)
cnt92 <- Read10X (data.dir="../3_clean_mtx/Ara_TR092_none/", gene.column=1)
obj91 <- CreateSeuratObject (counts=cnt91, min.features=1, min.cells=1)
obj92 <- CreateSeuratObject (counts=cnt92, min.features=1, min.cells=1)
obj91@meta.data$spl <- "ath_91"
obj92@meta.data$spl <- "ath_92"
obj <- merge (obj91, y=obj92, add.cell.ids=c("ath_91", "ath_92"))
obj <- JoinLayers (obj);
obj <- subset (obj, subset= nFeature_RNA>500 & nFeature_RNA<5000);

obj <- NormalizeData (obj);
obj <- FindVariableFeatures (obj);
obj <- ScaleData (obj);
obj <- RunPCA (obj);
obj <- FindNeighbors (obj, dims=1:20);
obj <- FindClusters (obj, resolution=0.5);
obj <- RunUMAP (obj, dims=1:20);
deg <- FindAllMarkers (Ara, only.pos=F);
saveRDS (obj, file="ath.gene_500-5000.resolution_0.5.rds");
write.table (deg, file="ath_91-92.all_deg.0.5.txt", row.names=T, col.names=T, sep='\t', quote=F);

pdf ("ath.gene_500-5000.resolution_0.5.pdf");
DimPlot (obj, group.by="spl");
DimPlot (obj, label=T);
dev.off ();
Ara <- readRDS("D:/sn-CUT&TAG-BGI_result/1_scRNA/1_scRNA/ath.gene_500-5000.resolution_0.5.anno.rds")
Araun <- subset(Ara, subset = cell_type != "unknown")

Idents(Ara) <- Ara@meta.data[["RNA_snn_res.0.5"]]

write.csv(
  deg,
  file = "E:/sn_CT/table/ara.degs.csv",
  row.names = TRUE,
  quote = FALSE
)