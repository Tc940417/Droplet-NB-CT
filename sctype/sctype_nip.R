lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)
source("D:/script/gene_sets_prepare.R")
source("D:/script/sctype_score.R")
nip.sctype <- readRDS(file="D:/sn-CUT&TAG-BGI_result/1_scRNA/1_scRNA/nip.gene_500-5000.resolution_0.25.rds")
scaled_matrix <- LayerData(nip.sctype, assay = "RNA", layer = "scale.data")

## load the customized database file
db_="D:/sn-CUT&TAG-BGI_result/1_scRNA/1_scRNA/nip.RNA.anno.marker_DB/ricemerge_DB -.xlsx"

## set the tissue type
tissue="Leaf"
gs_list=gene_sets_prepare(db_,tissue)
es.max = sctype_score(scRNAseqData = scaled_matrix, scaled = TRUE,
                      gs = gs_list$gs_positive, gs2 = NULL, gene_names_to_uppercase = FALSE)

cL_resutls = do.call("rbind", lapply(unique(nip.sctype@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(nip.sctype@meta.data[nip.sctype@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(nip.sctype@meta.data$seurat_clusters==cl)), 12)
}))

sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"


## check scores
print(sctype_scores[,1:3])


## add customclassif attributes
nip.sctype@meta.data$customclassif = ""
  for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,];
  nip.sctype@meta.data$customclassif[nip.sctype@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
  dim_label<-DimPlot(nip.sctype, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')}
print(dim_label)



saveRDS (nip.sctype, file="nip.sctype.anno.rds")
nip.sctype <- readRDS(file="D:/sn-CUT&TAG/1_scRNA/1_scRNA/nip.sctype.anno.rds")
  
