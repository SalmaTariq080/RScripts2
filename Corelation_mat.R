setwd("E:\\BreastCancer_Project")
library(tidyr)
data <- read.csv("exp_mat.csv", row.names = 1)
data2 <- read.csv("All_genes.csv")
Sig_genes <- read.csv("Sec_Tar.csv")
highlight_symbols <- Sig_genes$Symbol
tdata <- as.data.frame(t(data))
tdata2 <- as.data.frame(round(cor(tdata,use="pairwise.complete.obs"),2))
ESR1_data <- as.data.frame(tdata2$`ESR1|2099`)
ESR1_data$gene_id <- row.names(tdata2)
ESR1_data <- ESR1_data %>% separate(gene_id, c('Symbol', 'gene_id'))
ESR1_data$logfc <- data2$log2FoldChange
colnames(ESR1_data)[1] <- "ESR1_Cor"
CHRNA5_data <- as.data.frame(tdata2$`CHRNA5|1138`)
CHRNA5_data$gene_id <- row.names(tdata2)
CHRNA5_data <- CHRNA5_data %>% separate(gene_id, c('Symbol', 'gene_id'))
CHRNA5_data$logfc <- data2$log2FoldChange
colnames(CHRNA5_data)[1] <- "CHRNA5_Cor"
write.csv(data.frame(ID=row.names(tdata2),tdata2), "All_Cor.csv", quote = F, row.names = F)
write.csv(ESR1_data, "ESR1_Cor.csv", quote = F, row.names = F)
write.csv(CHRNA5_data, "CHRNA5_Cor.csv", quote = F, row.names = F)
pdf("Corfc_plot.pdf",width=10,height=10,paper='special')
plot(CHRNA5_data$logfc, CHRNA5_data$CHRNA5_Cor, main="Plot between correlation and foldchange", xlab="foldchange", ylab="CHRNA5_Cor")
highlight_points <- CHRNA5_data$Symbol %in% highlight_symbols
points(CHRNA5_data$logfc[highlight_points], CHRNA5_data$CHRNA5_Cor[highlight_points], col = "blue", pch = 16)
plot(ESR1_data$logfc, ESR1_data$ESR1_Cor, main="Plot between correlation and foldchange", xlab="foldchange", ylab="ESR1_Cor")
highlight_points <- ESR1_data$Symbol %in% highlight_symbols
points(ESR1_data$logfc[highlight_points], ESR1_data$ESR1_Cor[highlight_points], col = "blue", pch = 16)
dev.off()
library(ggplot2)
ggplot(ESR1_data, aes(x = logfc, y = ESR1_Cor, color = Symbol)) +
  geom_point() +  # Add points for the scatter plot
  labs(
    title = "Scatter Plot between Correlation and Fold Change",
    x = "Fold Change",
    y = "ESR1_Cor"
  ) +
  geom_point(data = subset(ESR1_data, Symbol %in% highlight_symbols), color = "red", size = 3) +
  guides(color = guide_legend(title = "Symbol"))  #
