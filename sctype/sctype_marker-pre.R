data <- read_excel("C:/Users/admin/Desktop/rice marker evi.xlsx")
result <- data %>%
  group_by(cell_type) %>%
  summarise(gene_id_combined = paste(gene_id, collapse = ","))
write.xlsx(result,file = "D:/sn-CUT&TAG/1_scRNA/1_scRNA/rice-evi_DB.xlsx")
marker <- read_excel("D:/sn-CUT&TAG/1_scRNA/1_scRNA/rice-evi_DB.xlsx")
