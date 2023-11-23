library(DESeq2)
library(IHW)
library(dplyr)
library(tidyr)
setwd("E:\\BreastCancer_Project")
data <- read.csv("final_countfile_pro.csv")
pheno_data <- read.csv("Phenotype_with_ER_status.csv")

#meta_data <- read.table("pancan_samples.txt", header=TRUE)
#deleting the first column as its just the numbering
pheno_data <- subset(pheno_data, select = -X )
#assigning the gene ids in the gene_id column as rownames and just deleting the first column so
# the table hasonly the count data 
#data <- data %>% separate(gene_id, c('Symbol', 'gene_id'))
rownames(data) <- data$gene_id
data <- subset(data, select = -gene_id )

#deleting rows where all entries are zero
data2 <- data[rowSums(data[])>0,]
#adding one to all entries in the data matrix and taking log2
data2 <- log2(data2+1)
#write.csv(data.frame(ID=row.names(data2),data2), "exp_mat.csv", quote = F, row.names = F)
#rm(data2)
rownames(pheno_data) <- pheno_data$patient_id
table(rownames(pheno_data) == pheno_data$patient_id)
all(rownames(pheno_data) %in% colnames(data))
all(rownames(pheno_data) == colnames(data))
data2 <- data2[, rownames(pheno_data)]
## DESeq2 Analysis
dds <- DESeqDataSetFromMatrix(countData = round(data2),
                              colData = pheno_data,
                              design = ~ ER.Status)
dds$ER.Status <- relevel(dds$ER.Status, ref = "Negative")
dds <- DESeq(dds)

resultsNames(dds)
#results use mean count filtering by default so i am specifying filtering 
#on the basis of independent hypothesis testing IHW package finds an optimal weighting 
#of the hypotheses that maximizes power while still controlling the FDR.
## DESeq2 results
res <- results(dds, filterFun = ihw, alpha = 0.05, name = "ER.Status_Positive_vs_Negative")
summary(res)
rm(data)
res_df <- as.data.frame(res)
#res_df$gene_id <- row.names(res_df)
#filtering the table on the basis of fold change and p value
res_Sigdf <- subset(res_df, (log2FoldChange < -1 | log2FoldChange > 1) & pvalue < 0.05)
res_Sig_dr <- subset(res_df, (log2FoldChange < -1) & pvalue < 0.05)
res_Sig_ur <- subset(res_df, ( log2FoldChange > 1) & pvalue < 0.05)

write.csv(data.frame(ID=row.names(res_Sigdf),res_Sigdf), "Sig_deg.csv", quote = F, row.names = F)
write.csv(data.frame(ID=row.names(res_df),res_df), "All_genes.csv", quote = F, row.names = F)
write.csv(data.frame(ID=row.names(res_Sig_dr),res_Sig_dr), "res_Sig_dr.csv", quote = F, row.names = F)
write.csv(data.frame(ID=row.names(res_Sig_ur),res_Sig_ur), "res_Sig_ur.csv", quote = F, row.names = F)
res_df <- read.csv("All_genes.csv")
library(EnhancedVolcano)
pdf("Volcano_plot_padj.pdf",width=10,height=10,paper='special')
EnhancedVolcano(res_df,
                lab = res_df$ID,
                x = 'log2FoldChange',
                y = 'pvalue', # change to padj for adjusted p value
                xlim = c(-5, 5),
                title = 'ER+ vs ER-',
                pCutoff = 0.05,
                FCcutoff = 1,
                #xlab = "Log2 fold change",  # Change the x-axis label
               # ylab = "Adjusted_p_value"
        )

dev.off()
#Creating PCA plot 
#Applying variance normalization: VST (Variance Stabilizing Transformation)
vst_dds <- vst(dds) 
pdf("PCA_plot.pdf",width=10,height=10,paper='special')
plotPCA(vst_dds, intgroup = c("ER.Status"))
dev.off()
rm(vst_dds)
