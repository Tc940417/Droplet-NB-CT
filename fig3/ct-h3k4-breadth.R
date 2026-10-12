
# ============================================================
# H3K4me3 Peak Breadth Analysis in Rice
# Association Between Cell-Type Marker Genes and Peak Breadth
# ============================================================

# Step 1. Load required libraries
library(Seurat)
library(Signac)
library(GenomicRanges)
library(dplyr)
library(ggplot2)

# Step 2. Set input and output directories
data_dir <- "data"
output_dir <- "results/peak_breadth"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Part 1. Identify Cell-Type Marker Genes from scRNA-seq
# ============================================================

# Step 3. Read annotated scRNA-seq data
RNA <- readRDS(
  file.path(data_dir, "nip_annotated.rds")
)

# Step 4. Retain control samples
ct_RNA <- subset(
  RNA,
  subset = group == "ctrl" & !is.na(cell_type) & cell_type != "Unknown"
)

# Join RNA layers for differential expression analysis
DefaultAssay(ct_RNA) <- "RNA"
ct_RNA <- JoinLayers(ct_RNA, assay = "RNA")

# Set cell type identities
Idents(ct_RNA) <- "cell_type"

# Step 5. Identify cell-type marker genes using MAST
deg <- FindAllMarkers(
  ct_RNA,
  assay = "RNA",
  test.use = "MAST",
  min.pct = 0.1,
  logfc.threshold = 0.5,
  latent.vars = "nCount_RNA",
  only.pos = TRUE
)

# Step 6. Filter significant marker genes
deg_filt <- deg %>%
  filter(
    avg_log2FC > 0.5,
    p_val_adj < 0.05
  ) %>%
  mutate(
    gene = gsub("-", "_", gene)
  )

# Save differential expression results
write.csv(
  deg,
  file = file.path(output_dir, "Rice_celltype_markers_MAST.csv"),
  row.names = FALSE
)

write.csv(
  deg_filt,
  file = file.path(output_dir, "Rice_celltype_markers_filtered.csv"),
  row.names = FALSE
)

# ============================================================
# Part 2. Calculate H3K4me3 Peak Breadth
# ============================================================

# Step 7. Read H3K4me3 Seurat object
h3k4 <- readRDS(
  file.path(data_dir, "ct.h3k4.rds")
)

DefaultAssay(h3k4) <- "ATAC"

# Step 8. Extract gene annotations
gene_anno <- Annotation(h3k4)
gene_anno <- gene_anno[gene_anno$type == "gene"]

# Step 9. Define promoter regions (+/- 2 kb around TSS)
promoter_regions <- promoters(
  gene_anno,
  upstream = 2000,
  downstream = 2000
)

# Step 10. Extract H3K4me3 peaks
peaks <- granges(h3k4)

peaks$breadth <- width(peaks)

# Step 11. Identify promoter-overlapping peaks
overlap <- findOverlaps(
  promoter_regions,
  peaks,
  ignore.strand = TRUE
)

# Step 12. Associate peaks with genes
peak_gene_df <- data.frame(
  gene = promoter_regions$gene_name[queryHits(overlap)],
  breadth = peaks$breadth[subjectHits(overlap)]
)

# Step 13. Calculate total peak breadth per gene
gene_breadth <- peak_gene_df %>%
  filter(!is.na(gene)) %>%
  group_by(gene) %>%
  summarise(
    peak_breadth = sum(breadth),
    .groups = "drop"
  )

# Remove zero or invalid values
gene_breadth <- gene_breadth %>%
  filter(is.finite(peak_breadth), peak_breadth > 0)

# Save gene-level peak breadth
write.csv(
  gene_breadth,
  file = file.path(output_dir, "Rice_H3K4me3_gene_peak_breadth.csv"),
  row.names = FALSE
)

# ============================================================
# Part 3. Compare All Genes and Marker Genes
# ============================================================

# Step 14. Prepare plot data
all_genes_breadth <- gene_breadth %>%
  mutate(group = "All genes")

deg_genes_breadth <- gene_breadth %>%
  filter(gene %in% deg_filt$gene) %>%
  mutate(group = "DEG genes")

plot_df <- bind_rows(
  all_genes_breadth,
  deg_genes_breadth
)

# Step 15. Plot peak breadth distributions
p1 <- ggplot(
  plot_df,
  aes(x = peak_breadth, fill = group)
) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 50,
    position = "identity",
    alpha = 0.55,
    color = "#3f3f3f",
    linewidth = 0.3
  ) +
  geom_density(
    aes(color = group),
    linewidth = 0.7,
    fill = NA
  ) +
  scale_x_log10() +
  scale_fill_manual(
    values = c(
      "All genes" = "#c97c6d",
      "DEG genes" = "#7f9d8c"
    )
  ) +
  scale_color_manual(
    values = c(
      "All genes" = "#9b5f54",
      "DEG genes" = "#5f7f6b"
    )
  ) +
  labs(
    x = "H3K4me3 peak breadth (bp)",
    y = "Density",
    fill = NULL,
    color = NULL
  ) +
  theme_classic(base_size = 12)

# Save distribution plot
ggsave(
  filename = file.path(
    output_dir,
    "Rice_H3K4me3_peak_breadth_distribution.pdf"
  ),
  plot = p1,
  device = "pdf",
  width = 8,
  height = 6
)

# ============================================================
# Part 4. Compare Peak Breadth Across Cell Types
# ============================================================

# Step 16. Associate marker genes with cell types
deg_anno <- deg_filt %>%
  select(gene, cluster) %>%
  distinct()

deg_genes_breadth <- gene_breadth %>%
  inner_join(deg_anno, by = "gene") %>%
  mutate(group = as.character(cluster))

all_genes_breadth <- gene_breadth %>%
  mutate(group = "All genes")

plot_df <- bind_rows(
  all_genes_breadth,
  deg_genes_breadth
)

# Step 17. Define cell type colors
cell_colors <- c(
  "Guard cell" = "#d69971",
  "Mesophyll" = "#64a776",
  "Epidermis" = "#b85292",
  "Large parenchyma" = "#537eb7",
  "Fiber" = "#efd2c9",
  "Vascular" = "#9f8d89",
  "All genes" = "#000000"
)

# Define plotting order
cell_type_order <- c(
  "Large parenchyma",
  "Mesophyll",
  "Epidermis",
  "Vascular",
  "Fiber",
  "Guard cell",
  "All genes"
)

plot_df$group <- factor(
  plot_df$group,
  levels = cell_type_order
)

# Step 18. Plot cumulative peak breadth distributions
p2 <- ggplot() +
  stat_ecdf(
    data = subset(plot_df, group != "All genes"),
    aes(x = peak_breadth, color = group),
    linewidth = 0.8
  ) +
  stat_ecdf(
    data = subset(plot_df, group == "All genes"),
    aes(x = peak_breadth, color = group),
    linewidth = 0.5,
    linetype = "dashed"
  ) +
  scale_x_continuous(
    breaks = c(0, 1000, 2000, 3000, 4000, 5000)
  ) +
  coord_cartesian(xlim = c(0, 5000)) +
  scale_color_manual(
    values = cell_colors
  ) +
  theme_classic(base_size = 12) +
  labs(
    x = "H3K4me3 peak breadth (bp)",
    y = "Cumulative fraction",
    color = "Cell type"
  )

# Step 19. Save cumulative distribution plot
ggsave(
  filename = file.path(
    output_dir,
    "Rice_H3K4me3_peak_breadth_ECDF.pdf"
  ),
  plot = p2,
  device = "pdf",
  width = 8,
  height = 6
)

# Step 20. Display plots
print(p1)
print(p2)
