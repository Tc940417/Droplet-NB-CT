
# ============================================================
# Rice Cell-Type Marker Gene Database Preparation
# Organize and Combine Marker Genes by Cell Type
# ============================================================

# Step 1. Load required libraries
library(readxl)
library(dplyr)
library(openxlsx)

# Step 2. Set input and output directories
data_dir <- "data"
output_dir <- "results"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Step 3. Read marker gene information
marker_data <- read_excel(
  file.path(data_dir, "rice marker evi.xlsx")
)

# Step 4. Combine marker genes by cell type
marker_result <- marker_data %>%
  group_by(cell_type) %>%
  summarise(
    gene_id_combined = paste(gene_id, collapse = ","),
    .groups = "drop"
  )

# Step 5. Save marker gene database
output_file <- file.path(output_dir, "rice-evi_DB.xlsx")

write.xlsx(
  marker_result,
  file = output_file,
  rowNames = FALSE,
  overwrite = TRUE
)

# Step 6. Read and display the generated database
marker <- read_excel(output_file)

print(marker)
