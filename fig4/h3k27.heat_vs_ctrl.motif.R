h3k27 <- readRDS("E:/sn_CT/RDS&h5ad/snCTaddRNAaddChromvar/h3k27.rds")
DefaultAssay(h3k27) <- 'chromvar'
h3k27$celltype <- paste(h3k27$batch, h3k27$cell_type, sep = "-")
Idents(h3k27) <- "celltype"



me.da <- FindMarkers(
  object = h3k27,
  ident.1 = 'heat_h3k27ac-Mesophyll',
  ident.2 = 'ctrl_h3k27ac-Mesophyll',
  only.pos = FALSE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)
lp.da <- FindMarkers(
  object = h3k27,
  ident.1 = 'heat_h3k27ac-Large parenchyma',
  ident.2 = 'ctrl_h3k27ac-Large parenchyma',
  only.pos = FALSE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)
ep.da <- FindMarkers(
  object = h3k27,
  ident.1 = 'heat_h3k27ac-Epidermis',
  ident.2 = 'ctrl_h3k27ac-Epidermis',
  only.pos = FALSE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)
vs.da <- FindMarkers(
  object = h3k27,
  ident.1 = 'heat_h3k27ac-Vascular',
  ident.2 = 'ctrl_h3k27ac-Vascular',
  only.pos = FALSE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)
fb.da <- FindMarkers(
  object = h3k27,
  ident.1 = 'heat_h3k27ac-Fiber',
  ident.2 = 'ctrl_h3k27ac-Fiber',
  only.pos = FALSE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)
gc.da <- FindMarkers(
  object = h3k27,
  ident.1 = 'heat_h3k27ac-Guard cell',
  ident.2 = 'ctrl_h3k27ac-Guard cell',
  only.pos = FALSE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

library(dplyr)
library(tibble)

da_list <- list(
  Mesophyll = me.da,
  `Large parenchyma` = lp.da,
  Epidermis = ep.da,
  Vascular = vs.da,
  Fiber = fb.da,
  `Guard cell` = gc.da
)

da_df <- lapply(names(da_list), function(ct) {
  da_list[[ct]] %>%
    rownames_to_column("feature") %>%
    mutate(
      cell_type = ct,
      neglog10_padj = -log10(p_val_adj),
      pct_diff = abs(pct.1 - pct.2)
    )
}) %>%
  bind_rows()



da_df <- da_df %>%
  mutate(
    padj_plot = ifelse(p_val_adj == 0, 1e-303, p_val_adj),
    neglog10_padj_plot = -log10(padj_plot)
  )






write.csv(
  me.da,
  file = "E:/sn_CT/table/H3K27ac_heat_vs_ctrl_Mesophyll_differential_motif_activity.csv",
  row.names = TRUE,
  quote = FALSE
)

write.csv(
  lp.da,
  file = "E:/sn_CT/table/H3K27ac_heat_vs_ctrl_Large_parenchyma_differential_motif_activity.csv",
  row.names = TRUE,
  quote = FALSE
)

write.csv(
  ep.da,
  file = "E:/sn_CT/table/H3K27ac_heat_vs_ctrl_Epidermis_differential_motif_activity.csv",
  row.names = TRUE,
  quote = FALSE
)

write.csv(
  vs.da,
  file = "E:/sn_CT/table/H3K27ac_heat_vs_ctrl_Vascular_differential_motif_activity.csv",
  row.names = TRUE,
  quote = FALSE
)

write.csv(
  fb.da,
  file = "E:/sn_CT/table/H3K27ac_heat_vs_ctrl_Fiber_differential_motif_activity.csv",
  row.names = TRUE,
  quote = FALSE
)

write.csv(
  gc.da,
  file = "E:/sn_CT/table/H3K27ac_heat_vs_ctrl_Guard_cell_differential_motif_activity.csv",
  row.names = TRUE,
  quote = FALSE
)
















