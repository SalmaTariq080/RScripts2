library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Proteome Profiling",
  data.type = "Protein Expression Quantification"
)
GDCdownload(query = query)
dat <- GDCprepare(query = query, save = TRUE, save.filename = "exp.rda")

colnames(P_Matrix)
