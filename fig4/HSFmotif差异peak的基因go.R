h3k4 <- readRDS("E:/sn_CT/RDS&h5ad/snCTaddRNAaddChromvar/h3k4.rds")
DefaultAssay(h3k4) <- "ATAC"
h3k4 <- RunTFIDF(h3k4)
h3k4$celltype <- paste(h3k4$batch, h3k4$cell_type, sep = "-")
Idents(h3k4) <- "celltype"
motif_mat <- GetMotifData(h3k4, slot = "data")
target_motifs <- c(
  "MA1664.3",
  "MA1758.2",
  "MA1759.2",
  "MA1760.1",
  "MA2019.2",
  "MA1667.3",
  "MA1665.3",
  "MA1666.3",
  "MA1761.2",
  "MA2020.2"
)
A.hsf_peaks <- rownames(motif_mat)[
  rowSums(motif_mat[, target_motifs, drop = FALSE]) > 0
]



hsf_da <- FindMarkers(
  h3k4,
  ident.1 = "heat_h3k4me3-Mesophyll",
  ident.2 = "ctrl_h3k4me3-Mesophyll",
  features = A.hsf_peaks,
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  only.pos = FALSE
)
hsf_da_sig <- hsf_da %>%
  dplyr::filter(
    avg_diff > 0.1,
    p_val_adj < 0.05
  )
sig_peaks <- rownames(hsf_da_sig)

peak2gene <- ClosestFeature(
  object = h3k4,
  regions = sig_peaks
)
sig_genes <- unique(peak2gene$gene_name)

up <- enrichGO(gene = sig_genes,
               OrgDb = org.Osativa.eg.db,
               keyType = "GID",    # 这里keyType视你的基因ID类型而定
               ont = "ALL",
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05)
up <- clusterProfiler::simplify(up,cutoff = 0.7,by = "p.adjust",select_fun = min,measure = "Wang")


dotplot(
  up,
  showCategory = 20
) +
  scale_color_gradientn(
    colors = c("#2166ac", "#67a9cf", "#fddbc7", "#b2182b")
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    legend.position = "right"
  )



go_df <- as.data.frame(up@result)
write.csv(
  go_df,
  file = "E:/sn_CT/table/ME_GO_enrichment_results.csv",
  row.names = FALSE,
  quote = FALSE
)
















