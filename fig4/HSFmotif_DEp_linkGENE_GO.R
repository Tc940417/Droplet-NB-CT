
# ============================================================
# HSF-Associated H3K4me3 Peak Analysis in Rice
# Differential Peak Analysis and GO Enrichment in Mesophyll
# ============================================================

# Step 1. Load required libraries
library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Osativa.eg.db)
library(enrichplot)

# Step 2. Set input and output directories
data_dir <- "data"
output_dir <- "results/HSF_GO"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Step 3. Read H3K4me3 Seurat object
h3k4 <- readRDS(
  file.path(data_dir, "h3k4.rds")
)

DefaultAssay(h3k4) <- "ATAC"

# Step 4. Normalize ATAC data using TF-IDF
h3k4 <- RunTFIDF(h3k4)

# Step 5. Define sample and cell type identities
h3k4$celltype <- paste(
  h3k4$batch,
  h3k4$cell_type,
  sep = "-"
)

Idents(h3k4) <- "celltype"

# ============================================================
# Part 1. Identify HSF Motif-Associated Peaks
# ============================================================

# Step 6. Define HSF motif IDs
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

# Step 7. Extract peak-by-motif occurrence matrix
motif_mat <- GetMotifData(
  h3k4,
  assay = "ATAC",
  slot = "data"
)

# Retain available HSF motifs
available_motifs <- intersect(
  target_motifs,
  colnames(motif_mat)
)

if (length(available_motifs) == 0) {
  stop("No target HSF motifs found in the motif matrix.")
}

if (length(available_motifs) < length(target_motifs)) {
  warning("Some target HSF motifs are missing.")
}

# Step 8. Identify peaks containing HSF motifs
A.hsf_peaks <- rownames(motif_mat)[
  Matrix::rowSums(
    motif_mat[, available_motifs, drop = FALSE]
  ) > 0
]

# Retain peaks present in the ATAC assay
A.hsf_peaks <- intersect(
  A.hsf_peaks,
  rownames(h3k4[["ATAC"]])
)

if (length(A.hsf_peaks) == 0) {
  stop("No HSF-associated peaks found.")
}

# ============================================================
# Part 2. Differential H3K4me3 Peak Analysis
# ============================================================

# Step 9. Compare heat-treated and control Mesophyll cells
hsf_da <- FindMarkers(
  object = h3k4,
  ident.1 = "heat_h3k4me3-Mesophyll",
  ident.2 = "ctrl_h3k4me3-Mesophyll",
  assay = "ATAC",
  slot = "data",
  features = A.hsf_peaks,
  test.use = "wilcox",
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  logfc.threshold = 0,
  min.pct = 0,
  only.pos = FALSE
)

# Step 10. Filter significantly increased peaks
hsf_da_sig <- hsf_da %>%
  filter(
    avg_diff > 0.1,
    p_val_adj < 0.05
  )

sig_peaks <- rownames(hsf_da_sig)

# Step 11. Save differential peak results
write.csv(
  hsf_da,
  file = file.path(output_dir, "HSF_differential_peaks.csv"),
  row.names = TRUE,
  quote = FALSE
)

write.csv(
  hsf_da_sig,
  file = file.path(output_dir, "HSF_significant_peaks.csv"),
  row.names = TRUE,
  quote = FALSE
)

if (length(sig_peaks) == 0) {
  stop("No significant HSF-associated peaks identified.")
}

# ============================================================
# Part 3. Peak-to-Gene Annotation
# ============================================================

# Step 12. Annotate significant peaks with nearest genes
peak2gene <- ClosestFeature(
  object = h3k4,
  regions = sig_peaks
)

# Step 13. Extract associated genes
sig_genes <- unique(
  na.omit(peak2gene$gene_name)
)

# Step 14. Save peak-to-gene annotations
write.csv(
  peak2gene,
  file = file.path(output_dir, "HSF_peak_to_gene.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  data.frame(gene = sig_genes),
  file = file.path(output_dir, "HSF_associated_genes.csv"),
  row.names = FALSE,
  quote = FALSE
)

# ============================================================
# Part 4. GO Enrichment Analysis
# ============================================================

# Step 15. Check gene identifier compatibility
key_type <- "GID"

if (!key_type %in% AnnotationDbi::keytypes(org.Osativa.eg.db)) {
  stop("GID is not supported by org.Osativa.eg.db.")
}

# Retain genes recognized by the annotation database
valid_keys <- AnnotationDbi::keys(
  org.Osativa.eg.db,
  keytype = key_type
)

go_genes <- intersect(sig_genes, valid_keys)

message("Associated genes: ", length(sig_genes))
message("GO-matched genes: ", length(go_genes))

if (length(go_genes) == 0) {
  stop("No gene IDs match the GO annotation database.")
}

# Step 16. Perform GO enrichment
up <- enrichGO(
  gene = go_genes,
  OrgDb = org.Osativa.eg.db,
  keyType = key_type,
  ont = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = FALSE
)

# Step 17. Save GO enrichment results
go_df <- as.data.frame(up)

write.csv(
  go_df,
  file = file.path(output_dir, "HSF_GO_enrichment_results.csv"),
  row.names = FALSE,
  quote = FALSE
)

# ============================================================
# Part 5. GO Enrichment Visualization
# ============================================================

# Step 18. Generate GO enrichment dotplot
if (nrow(go_df) > 0) {

  p <- enrichplot::dotplot(
    up,
    showCategory = 20
  ) +
    scale_color_gradientn(
      colours = c(
        "#2166ac",
        "#67a9cf",
        "#fddbc7",
        "#b2182b"
      )
    ) +
    theme_classic() +
    theme(
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      legend.position = "right"
    )

  # Step 19. Save GO enrichment plot
  ggsave(
    filename = file.path(
      output_dir,
      "HSF_GO_enrichment_dotplot.pdf"
    ),
    plot = p,
    device = "pdf",
    width = 10,
    height = 8
  )

  print(p)

} else {
  warning("No significant GO terms were identified.")
}
