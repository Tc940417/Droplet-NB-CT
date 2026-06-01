h3k27 <- readRDS("E:/sn_CT/RDS&h5ad/snCTaddRNAaddChromvar/ct.h3k27.rds")
motif_avg <- AverageExpression(
  h3k27,
  assays = "chromvar",
  slot = "data",
  group.by = "cell_type"
)
motif_avg_matrix <- motif_avg$chromvar
motif_avg_matrix <- as.matrix(motif_avg$chromvar)

library(ComplexHeatmap)
library(circlize)

Heatmap(
  motif_avg_matrix,
  name = "Deviation",
  col = colorRamp2(c(-2, 0, 2), c("#000A3D", "white", "#D7301F")),
  show_row_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE
)
top10_list <- lapply(
  colnames(motif_avg_matrix),
  function(ct) {
    vals <- motif_avg_matrix[, ct]
    top_motifs <- names(sort(vals, decreasing = TRUE))[1:10]
    data.frame(
      cell_type = ct,
      motif = top_motifs,
      deviation = vals[top_motifs]
    )
  }
)

top10_motif_df <- do.call(rbind, top10_list)
top10_motif_df

df_split <- split(top10_motif_df, top10_motif_df$cell_type)

MotifPlot(
  object = h3k27,
  motifs = df_split$motif,
  assay = 'ATAC'
)



write.csv(
  motif_avg_matrix,
  file = "E:/sn_CT/table/h3k27motif_avg_matrix.csv",
  row.names = TRUE,
  quote = FALSE
)





