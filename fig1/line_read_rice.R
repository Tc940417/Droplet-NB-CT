df3 = read.table("D:/sn-CT_broadpeaks/3/singlecell.csv",header=TRUE,sep=',')
rds3 = snCT_3

df4 = read.table("D:/sn-CT_broadpeaks/4/singlecell.csv",header=TRUE,sep=',')
rds4 = snCT_4

df5 = read.table("D:/sn-CT_broadpeaks/5/singlecell.csv",header=TRUE,sep=',')
rds5 = snCT_5

df6 = read.table("D:/sn-CT_broadpeaks/6/singlecell.csv",header=TRUE,sep=',')
rds6 = snCT_6

df7 = read.table("D:/sn-CT_broadpeaks/7/singlecell.csv",header=TRUE,sep=',')
rds7 = snCT_7

df8 = read.table("D:/sn-CT_broadpeaks/8/singlecell.csv",header=TRUE,sep=',')
rds8 = snCT_8


df_1 = subset(df8, select = c(Cell, totalFrags))
df_2 = subset(df_1, df_1$Cell %in%colnames(rds8))
df_3 <- df_2[match(colnames(rds8), df_2$Cell), ]
rds8$nReads <- df_3$totalFrags

mean(rds3@meta.data$nFeature_CT[which(rds3@meta.data$nReads < 410000 & rds3@meta.data$nReads > 0)])

###nFeature_CT###
# 初始化结果向量
means <- c()

# 构建不规则递增的上限区间：
# 6000–10000: 每1000
# 10000–100000: 每10000
# 100000–400000: 每50000
upper_limits <- c(
  seq(7000, 10000, by = 1000),
  seq(11000, 100000, by = 10000),
  seq(150000, 450000, by = 20000)
)

# 循环计算每个区间的 mean(nFeature_CT)
for (upper in upper_limits) {
  # 提取 nReads 在 (0, upper] 区间的细胞
  cells_in_range <- which(rds8@meta.data$nReads > 0 & rds8@meta.data$nReads <= upper)
  
  # 计算均值（如无数据则为 NA）
  if (length(cells_in_range) > 0) {
    mean_value <- mean(rds8@meta.data$nFeature_CT[cells_in_range])
  } else {
    mean_value <- NA
  }
  
  # 存储结果
  means <- c(means, mean_value)
}

# 添加区间标签
names(means) <- paste0("(0, ", upper_limits, "]")

# 转换为 data.frame
df <- data.frame(
  Range = names(means),
  Mean_nFeature_CT = means
)

# 保存为 CSV 文件（当前工作目录下）
write.csv(df, file = "means_nFeature_CT8_by_nReads_custom.csv", row.names = FALSE)



###nCount_CT###
# 初始化结果向量
means <- c()

# 构建不规则递增的上限区间：
# 6000–10000: 每1000
# 10000–100000: 每10000
# 100000–400000: 每50000
upper_limits <- c(
  seq(7000, 10000, by = 1000),
  seq(11000, 100000, by = 10000),
  seq(150000, 450000, by = 20000)
)

# 循环计算每个区间的 mean(nCount_CT)
for (upper in upper_limits) {
  # 提取 nReads 在 (0, upper] 区间的细胞
  cells_in_range <- which(rds8@meta.data$nReads > 0 & rds8@meta.data$nReads <= upper)
  
  # 计算均值（如无数据则为 NA）
  if (length(cells_in_range) > 0) {
    mean_value <- mean(rds8@meta.data$nCount_CT[cells_in_range])
  } else {
    mean_value <- NA
  }
  
  # 存储结果
  means <- c(means, mean_value)
}

# 添加区间标签
names(means) <- paste0("(0, ", upper_limits, "]")

# 转换为 data.frame
df <- data.frame(
  Range = names(means),
  Mean_nCount_CT = means
)

# 保存为 CSV 文件（当前工作目录下）
write.csv(df, file = "means_nCount_CT8_by_nReads_custom.csv", row.names = FALSE)






































