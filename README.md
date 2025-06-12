# muggle_genomics
#TPM_retrieval_TCGA:ForWT_p53_cohort
library(TCGAbiolinks)
library(SummarizedExperiment)
data <- read.csv("p53WT_cohort_Barcodes.csv")
TP_barcode <- data$Barcodes
q1 <- TCGAquery_SampleTypes(TP_barcode,"TP") 
TP_query <- GDCquery(
  project = c("TCGA-COAD"),
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts",
  barcode=TP_barcode)
GDCdownload(TP_query) 
data <- GDCprepare(TP_query)
TPM <- assay(data,4) %>% data.frame()
head(TPM)
dim(TPM)
write.csv(TPM,"WT_TPM_counts.csv")
