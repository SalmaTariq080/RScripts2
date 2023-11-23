library(dplyr)
# reading the data and removing all na values
if (any(is.na(data))) {
  # Remove missing values
  P_Matrix <- na.omit(data)
}
###########################################
#i omitted na values but so much values are completely remove if the percentage of na values is 0.001088 it still got removed

#################################################################
P_Matrix <- data
rownames(P_Matrix) <- P_Matrix[,1]
P_Matrix <- select(P_Matrix, -c(1:5))
###################################################################
# so tried the following code to omit NA values where percentage of NA values is graeter than 80%
# No values got deleted as none of the entry have NA values greater than 80%
# Calculate the percentage of zeros in each row
row_zero_percent <- rowMeans(data == 0, na.rm = TRUE)

# Keep only rows where less than 80% of the entries are zero
threshold <- 0.8
filtered_data <- data[row_zero_percent < threshold, ]
rm(filtered_data)
###########################################################3
clin_data <- read.csv("Clinical_data.csv")
#removing the last bit from the column names of P_Matrix
colnames(P_Matrix) <- gsub("-01A$|-01B$|-11A$|-11B$|-06A$", "", colnames(P_Matrix))
#getting the column names of p_matrix
colnames_P_Matrix <- colnames(P_Matrix)
# getting column names of pamtrix if they are also present in clin_data
colnames_P_Matrix <- colnames_P_Matrix[colnames_P_Matrix %in% clin_data$bcr_patient_barcode]
# getting the subset of P_Matrix 
df2_filtered <- P_Matrix[, colnames_P_Matrix]
clin_data_filtered <- clin_data[clin_data$bcr_patient_barcode %in% colnames(df2_filtered), ]
clin_data_filtered <- clin_data_filtered[clin_data_filtered$er_status_by_ihc != "[Not Evaluated]", ]
clin_data_filtered <- clin_data_filtered[clin_data_filtered$er_status_by_ihc != "Indeterminate", ]
rownames(clin_data_filtered) <- clin_data_filtered[,1]
df2_filtered_filtered <- df2_filtered %>%
  select(cols = all_of(clin_data_filtered$bcr_patient_barcode))

df2_filtered_filtered <- df2_filtered %>%
  select(matches(paste0("^", clin_data_filtered$bcr_patient_barcode, "$")))
write.csv(df2_filtered_filtered, "RPP_Exp_Mat.csv", quote = F, row.names = T)
write.csv(clin_data_filtered, "RPP_Clin_Data.csv", quote = F, row.names = F)
