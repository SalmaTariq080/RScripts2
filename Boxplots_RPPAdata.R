setwd("E:\\BreastCancer_Project")
clin_data <- read.csv("RPP_Clin_Data.csv")
P_Matrix <- read.csv("RPP_Exp_Mat.csv")
rownames(P_Matrix) <- P_Matrix[,1]
P_Matrix <- select(P_Matrix, -c(1))
ERAlpha <- "AGID00335"
ERRAlpha <- "AGID00405"
ERAlphas118 <- "AGID00029" 
ERAlpha_data <- P_Matrix[ERAlpha, , drop = FALSE]
ERRAlpha_data <- P_Matrix[ERRAlpha, , drop = FALSE]
ERAlphas118_data <- P_Matrix[ERAlphas118, , drop = FALSE]
combined_data <- rbind(ERAlpha_data, ERRAlpha_data, ERAlphas118_data)
t_combined_data <- as.data.frame(t(combined_data))
bcr_patient_barcode <- rownames(t_combined_data)
t_combined_data <- cbind(bcr_patient_barcode,t_combined_data)
clin_data$bcr_patient_barcode <- gsub("-", ".", clin_data$bcr_patient_barcode)
merged_data <- merge(t_combined_data, clin_data, by = "bcr_patient_barcode")
#########################################################
library(ggplot2)
library(RColorBrewer)
my_palette <- brewer.pal(2, "Set3")
pdf("Protein_boxplots.pdf",width=20,height=10,paper='special')
ggplot(merged_data, aes(x = er_status_by_ihc, y = AGID00335, fill = er_status_by_ihc)) +
  geom_boxplot() +
  scale_fill_manual(values = my_palette) +
  labs(title = "Boxplot for ERAlpha Expression in ER+ and ER- ", y = "Expression Level") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(merged_data, aes(x = er_status_by_ihc, y = AGID00405, fill = er_status_by_ihc)) +
  geom_boxplot() +
  scale_fill_manual(values = my_palette) +
  labs(title = "Boxplot for ERRAlpha Expression in ER+ and ER- ", y = "Expression Level") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(merged_data, aes(x = er_status_by_ihc, y = AGID00029, fill = er_status_by_ihc)) +
  geom_boxplot() +
  scale_fill_manual(values = my_palette) +
  labs(title = "Boxplot for ERAlphas118 Expression in ER+ and ER- ", y = "Expression Level") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()