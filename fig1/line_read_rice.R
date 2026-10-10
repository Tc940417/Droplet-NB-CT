library(Seurat)

# Load cell metadata and the single-cell histone modification Seurat object.
metadata <- read.csv(
  "path/to/singlecell.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)

rds <- readRDS("path/to/snCT.rds")

# Match fragment counts to the cell order in the Seurat object.
cell_index <- match(colnames(rds), metadata$Cell)
rds$nReads <- metadata$totalFrags[cell_index]

# Mean detected peaks for cells with 0 < nReads < 410000.
cells <- which(rds$nReads > 0 & rds$nReads < 410000)
mean(rds$nFeature_CT[cells])

# Preserve the upper limits in the original code.
upper_limits <- c(
  seq(7000, 10000, by = 1000),
  seq(11000, 100000, by = 10000),
  seq(150000, 450000, by = 20000)
)

# Calculate cumulative means for cells with 0 < nReads <= upper.
calculate_means <- function(metric) {
  means <- vapply(
    upper_limits,
    function(upper) {
      cells <- which(rds$nReads > 0 & rds$nReads <= upper)

      if (length(cells) == 0L) {
        return(NA_real_)
      }

      mean(rds[[metric]][cells, 1])
    },
    numeric(1)
  )

  result <- data.frame(
    Range = paste0("(0, ", upper_limits, "]"),
    Mean = means
  )
  colnames(result)[2] <- paste0("Mean_", metric)

  result
}

write.csv(
  calculate_means("nFeature_CT"),
  file = "means_nFeature_CT_by_nReads_custom.csv",
  row.names = FALSE
)

write.csv(
  calculate_means("nCount_CT"),
  file = "means_nCount_CT_by_nReads_custom.csv",
  row.names = FALSE
)
