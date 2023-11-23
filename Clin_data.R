library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
query <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Clinical",
  data.type = "Clinical Supplement", 
  data.format = "BCR Biotab"
)
GDCdownload(query)
clinical_tab_all <- GDCprepare(query)
names(clinical_tab_all)
Clin_pat <- dplyr::glimpse(clinical_tab_all$clinical_patient_brca)
subset_data <- Clin_pat[, c("bcr_patient_barcode", "er_status_by_ihc")]
write.csv(data.frame(subset_data), "Clinical_data.csv", quote = F, row.names = F)
