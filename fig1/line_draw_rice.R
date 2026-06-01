library(ggplot2)
library(reshape2)



count <- read.csv("means_nCount_by_nReads_custom.csv")
count_long <- melt(count, id.vars = "X", variable.name = "Group", value.name = "Value")

custom_colors <- c("#df5734", "#e97e61", "#64a776", "#8fc19a")
pdf(paste("fragments-reads_rice.pdf",sep=""),width = 8, height = 5)
ggplot(count_long, aes(x = X, y = Value, color = Group)) +
  geom_line(size = 1.1) + geom_point(size = 2) + scale_color_manual(values = custom_colors) +
  labs(title = "Mean of detected fragments", x = "Number of reads", y = "Mean of detected fragments", color = "Group") + theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,color = 'black',face = "bold", size = 12),
        axis.text.y = element_text(color = 'black', face = "bold", size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color="black",size = 1.2, linetype="solid"),
        panel.spacing = unit(0.12, "cm"),
        plot.title = element_text(hjust = 0.5,size = 14, face = "bold") 
  )
dev.off()



peak <- read.csv("means_nFeature_by_nReads_custom.csv")
peak_long <- melt(peak, id.vars = "X", variable.name = "Group", value.name = "Value")

custom_colors <- c("#df5734", "#e97e61", "#64a776", "#8fc19a")
pdf(paste("peaks-reads_rice.pdf",sep=""),width = 8, height = 5)
ggplot(peak_long, aes(x = X, y = Value, color = Group)) +
  geom_line(size = 1.1) + geom_point(size = 2) + scale_color_manual(values = custom_colors) +
  labs(title = "Mean of detected peaks", x = "Number of reads", y = "Mean of detected peaks", color = "Group") + theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,color = 'black',face = "bold", size = 12),
        axis.text.y = element_text(color = 'black', face = "bold", size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color="black",size = 1.2, linetype="solid"),
        panel.spacing = unit(0.12, "cm"),
        plot.title = element_text(hjust = 0.5,size = 14, face = "bold") 
  )
dev.off()



