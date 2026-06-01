library(dplyr)
library(Signac)
library(Seurat)
rna.me <- read.csv("E:/sn_CT/table/rice_snRNA_DEG/Mesophyll_case_vs_ctrl.csv",row.names = 1)
rna.lp <- read.csv("E:/sn_CT/table/rice_snRNA_DEG/LargeParenchyma_case_vs_ctrl.csv",row.names = 1)
rna.gc <- read.csv("E:/sn_CT/table/rice_snRNA_DEG/GuardCell_case_vs_ctrl.csv",row.names = 1)
rna.vs <- read.csv("E:/sn_CT/table/rice_snRNA_DEG/Vascular_case_vs_ctrl.csv",row.names = 1)
rna.ep <- read.csv("E:/sn_CT/table/rice_snRNA_DEG/Epidermis_case_vs_ctrl.csv",row.names = 1)
rna.fb <- read.csv("E:/sn_CT/table/rice_snRNA_DEG/Fiber_case_vs_ctrl.csv",row.names = 1)
DefaultAssay(h3k4)<- "ATAC"
DefaultAssay(h3k27)<- "ATAC"


celltypes <- c("me","lp","gc","vs","ep","fb")

rna_list <- list(me=rna.me, lp=rna.lp, gc=rna.gc, vs=rna.vs, ep=rna.ep, fb=rna.fb)
peak_list <- list(me=k27.me,  lp=k27.lp,  gc=k27.gc,  vs=k27.vs,  ep=k27.ep,  fb=k27.fb)

rna_f_list <- lapply(rna_list, function(df){
  df %>%
    filter(p_val_adj <0.01, abs(avg_log2FC) > 1.5) %>%
    rownames_to_column("gene")
})

rna_step1_summary <- data.frame(
  celltype = names(rna_f_list),
  n_deg = sapply(rna_f_list, nrow)
)

rna_step1_summary


peak_f_list <- lapply(peak_list, function(df){
  df %>%
    filter(p_val< 0.005, abs(avg_log2FC) > 1) %>%
    rownames_to_column("peak")
})

peak_step2_summary <- data.frame(
  celltype = names(peak_f_list),
  n_diffpeak = sapply(peak_f_list, nrow)
)

peak_step2_summary

#######################################################
peak_anno_list <- lapply(peak_f_list, function(df){
  
  peak.gr <- StringToGRanges(df$peak)
  
  cf <- ClosestFeature(
    object = h3k27,
    regions = peak.gr
  )
  
  # cf 里一定包含：
  # query_region, gene_name, distance
  
  df$gene <- cf$gene_name[
    match(df$peak, cf$query_region)
  ]
  
  df$distance <- cf$distance[
    match(df$peak, cf$query_region)
  ]
  
  df
})
peak_anno_list <- lapply(peak_anno_list, function(df){
  df %>%
    mutate(
      gene = gsub("_", "-", gene)
    )
})
peak_step3_summary <- data.frame(
  celltype = names(peak_anno_list),
  n_peak_total = sapply(peak_anno_list, nrow),
  n_peak_with_gene = sapply(peak_anno_list, function(df) sum(!is.na(df$gene)))
)

peak_step3_summary

#############################################################
peak_gene_list <- lapply(peak_anno_list, function(df){
  df %>%
    filter(!is.na(gene)) %>%
    group_by(gene) %>%
    summarise(
      peak_logFC = mean(avg_log2FC),
      n_peak = n(),
      mean_distance = mean(distance),
      .groups = "drop"
    ) %>%
    filter(abs(peak_logFC) > 1)
})

peak_step4_summary <- data.frame(
  celltype = names(peak_gene_list),
  n_gene_from_peaks = sapply(peak_gene_list, nrow)
)

peak_step4_summary

###############################################################

merged_list <- mapply(function(rna_df, peak_gene_df){
  inner_join(
    rna_df %>% select(gene, rna_logFC = avg_log2FC),
    peak_gene_df %>% select(gene, peak_logFC),
    by = "gene"
  )
}, rna_f_list, peak_gene_list, SIMPLIFY = FALSE)


merge_step5_summary <- data.frame(
  celltype = names(merged_list),
  n_overlap_gene = sapply(merged_list, nrow)
)

merge_step5_summary

#################################################################


cor_results <- do.call(rbind, lapply(names(merged_list), function(ct){
  df <- merged_list[[ct]]
  
  if(nrow(df) < 10){
    return(data.frame(celltype=ct, rho=NA, pvalue=NA, n_gene=nrow(df)))
  }
  
  res <- cor.test(df$rna_logFC, df$peak_logFC, method="spearman")
  
  data.frame(
    celltype = ct,
    rho = as.numeric(res$estimate),
    pvalue = res$p.value,
    n_gene = nrow(df)
  )
}))

cor_results %>% arrange(desc(abs(rho)))


saveRDS(merged_list,"E:/sn_CT/table/snCT_celltype_marker/snct&snRNA-cor/h3k27.rds")





library(ggplot2)

ct <- "lp"
df <- merged_list[[ct]]

# 计算相关系数（用于标注）
cor_val <- cor(df$peak_logFC, df$rna_logFC, method = "spearman")

ggplot(df, aes(x = peak_logFC, y = rna_logFC)) +
  geom_point(size = 2, alpha = 0.7, color = "#537eb7") +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black",
    fill = "black",
    alpha = 0.25
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "diff H3K4me3 (avg_log2FC)",
    y = "RNA DEGs (avg_log2FC)",
    title = paste0(ct, "  Correlation: ", round(cor_val, 2))
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70")



ct <- "me"
df <- merged_list[[ct]]

# 计算相关系数（用于标注）
cor_val <- cor(df$peak_logFC, df$rna_logFC, method = "spearman")

ggplot(df, aes(x = peak_logFC, y = rna_logFC)) +
  geom_point(size = 2, alpha = 0.7, color = "#64a776") +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black",
    fill = "black",
    alpha = 0.25
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "diff H3K4me3 (avg_log2FC)",
    y = "RNA DEGs (avg_log2FC)",
    title = paste0(ct, "  Correlation: ", round(cor_val, 2))
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70")





ct <- "ep"
df <- merged_list[[ct]]

# 计算相关系数（用于标注）
cor_val <- cor(df$peak_logFC, df$rna_logFC, method = "spearman")

ggplot(df, aes(x = peak_logFC, y = rna_logFC)) +
  geom_point(size = 2, alpha = 0.7, color = "#b85292") +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black",
    fill = "black",
    alpha = 0.25
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "diff H3K4me3 (avg_log2FC)",
    y = "RNA DEGs (avg_log2FC)",
    title = paste0(ct, "  Correlation: ", round(cor_val, 2))
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70")



ct <- "gc"
df <- merged_list[[ct]]

# 计算相关系数（用于标注）
cor_val <- cor(df$peak_logFC, df$rna_logFC, method = "spearman")

ggplot(df, aes(x = peak_logFC, y = rna_logFC)) +
  geom_point(size = 2, alpha = 0.7, color = "#d69971") +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black",
    fill = "black",
    alpha = 0.25
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "diff H3K4me3 (avg_log2FC)",
    y = "RNA DEGs (avg_log2FC)",
    title = paste0(ct, "  Correlation: ", round(cor_val, 2))
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70")





ct <- "fb"
df <- merged_list[[ct]]

# 计算相关系数（用于标注）
cor_val <- cor(df$peak_logFC, df$rna_logFC, method = "spearman")

ggplot(df, aes(x = peak_logFC, y = rna_logFC)) +
  geom_point(size = 2, alpha = 0.7, color = "#efd2c9") +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black",
    fill = "black",
    alpha = 0.25
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "diff H3K4me3 (avg_log2FC)",
    y = "RNA DEGs (avg_log2FC)",
    title = paste0(ct, "  Correlation: ", round(cor_val, 2))
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70")




ct <- "vs"
df <- merged_list[[ct]]

# 计算相关系数（用于标注）
cor_val <- cor(df$peak_logFC, df$rna_logFC, method = "spearman")

ggplot(df, aes(x = peak_logFC, y = rna_logFC)) +
  geom_point(size = 2, alpha = 0.7, color = "#9f8d89") +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black",
    fill = "black",
    alpha = 0.25
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "diff H3K4me3 (avg_log2FC)",
    y = "RNA DEGs (avg_log2FC)",
    title = paste0(ct, "  Correlation: ", round(cor_val, 2))
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70")

