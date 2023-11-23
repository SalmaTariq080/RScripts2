setwd("E:\\BreastCancer_Project")
library(tidyr)
library(biomaRt)
Sig_genes <- read.csv("Sec_Tar.csv")
gene_symbols <- Sig_genes$Symbol

mart=useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
mapping <- getBM(attributes=c("hgnc_id","uniprot_gn_symbol"), filters = "uniprot_gn_symbol", mart=mart, values=gene_symbols, uniqueRows=TRUE, bmHeader = T)
mapping <- mapping %>% separate("HGNC ID", c('HGNC', 'HGNC_ID'))
mapping <- subset(mapping, select = -HGNC )
colnames(mapping)[2] ="UniProt_Symbol"
data2 <- mapping[mapping$UniProt_Symbol %in% unique(Sig_genes$Symbol),]

unique_genes_in_sig <- setdiff(unique(Sig_genes$Symbol), unique(mapping$UniProt_Symbol))
mapping2 <- getBM(attributes=c("hgnc_id","uniprot_gn_symbol"), filters = "go", mart=mart, values=unique_genes_in_sig, uniqueRows=TRUE, bmHeader = T)
#################################################################################
# Load libraries
library(biomaRt)
library(dplyr)
library(tidyr)

# Function to retrieve gene information from BioMart
retrieve_gene_info <- function(gene_symbols) {
  mart <- useMart("ensembl")
  mart <- useDataset("hsapiens_gene_ensembl", mart)
  
  mapping <- getBM(
    attributes = c("hgnc_id", "uniprot_gn_symbol"),
    filters = "hgnc_symbol",
    values = gene_symbols,
    mart = mart,
    uniqueRows = TRUE,
    bmHeader = TRUE
  )
  
  # Process and return the mapping data
  mapping <- mapping %>% separate("HGNC ID", c('HGNC', 'HGNC_ID'))
  return(mapping)
}

# Main code
setwd("E:\\BreastCancer_Project")
Sig_genes <- read.csv("Sec_Tar.csv")
gene_symbols <- Sig_genes$Symbol

# Retrieve gene information
mapping <- retrieve_gene_info(gene_symbols)

# Handle unique genes
unique_genes_in_sig <- setdiff(unique(Sig_genes$Symbol), unique(mapping$hgnc_symbol))

if (length(unique_genes_in_sig) > 0) {
  mapping2 <- retrieve_gene_info(unique_genes_in_sig)
} else {
  mapping2 <- NULL
}
