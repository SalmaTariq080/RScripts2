library(DESeq2)
library(tidyverse)
library(dplyr)
library(tidyr)
library(ggpubr)
setwd("E:\\BreastCancer_Project")
clin_data <- read.csv("RPP_Clin_Data.csv")
P_Matrix <- read.csv("RPP_Exp_Mat.csv")
rownames(P_Matrix) <- P_Matrix[,1]
P_Matrix <- select(P_Matrix, -c(1))
t_P_Matrix <- as.data.frame(t(P_Matrix))
ERAlpha <- "AGID00335"
ERRAlpha <- "AGID00405"
ERAlphas118 <- "AGID00029" 
ERAlpha_data <- P_Matrix[ERAlpha, , drop = FALSE]
ERRAlpha_data <- P_Matrix[ERRAlpha, , drop = FALSE]
ERAlphas118_data <- P_Matrix[ERALPHApS118, , drop = FALSE]
#t_P_Matrix$bcr_patient_barcode <- row.names(t_P_Matrix)
###########################################################################
mrna_matrix <- read.csv("exp_mat.csv")
rownames(mrna_matrix) <- mrna_matrix[,1]
mrna_matrix <- select(mrna_matrix, -c(1))
ESR1 <- "ESR1|2099"
ESR1_data <- mrna_matrix[ESR1, , drop = FALSE]
common_genes <- intersect(names(ESR1_data), names(ERAlpha_data))
common_ESR1_data <- ESR1_data[, common_genes]
common_ERAlpha_data <- ERAlpha_data[, common_genes]
common_ERRAlpha_data <- ERRAlpha_data[, common_genes]
common_ERAlphas118_data <- ERAlphas118_data[, common_genes]
combined_data <- rbind(common_ESR1_data, common_ERAlpha_data, common_ERRAlpha_data, common_ERAlphas118_data)
t_combined_data <- as.data.frame(t(combined_data))
z_scores <- as.data.frame(scale(t_combined_data))
colnames(z_scores) <- c("ESR1","ERAlpha","ERRAlpha", "ERAlphapS118" )
pdf("Protein_gene_cor.pdf",width=20,height=10,paper='special')
ggscatter(z_scores, x = "ESR1", y = "ERAlpha", 
          add = "reg.line", conf.int = FALSE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "ESR1", ylab = "ERAlpha", color = "blue",
          fill = "gray", shape = 19,
          size = 2)
ggscatter(z_scores, x = "ESR1", y = "ERRAlpha", 
          add = "reg.line", conf.int = FALSE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "ESR1", ylab = "ERRAlpha", color = "blue",
          fill = "gray", shape = 19,
          size = 2)
ggscatter(z_scores, x = "ESR1", y = "ERAlphapS118", 
          add = "reg.line", conf.int = FALSE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "ESR1", ylab = "ERAlphapS118", color = "blue",
          fill = "gray", shape = 19,
          size = 2)
dev.off()
# AGID00405 for ERRAlpha
# AGID00335 for ERALPHA
# AGID00029 for ERALPHA_pS118
# ESR1|2099


