#!/usr/bin/env R

# Author: Sean Maden
#
# Read in TPM expression data. Data sources:
# 
# * "GSE107011_Processed_data_TPM.txt", GEO record.
#
#
#
#

libv <- c("SummarizedExperiment", "biomaRt")
sapply(libv, library, character.only = TRUE)

#----------------------------
# read and process expression
#----------------------------

# read geo expression into SummarizedExperiment
geoTpmTable <- "./data/monaco_et_al_2019/geo/GSE107011_Processed_data_TPM.txt"
tpm <- read.table(geoTpmTable)
newSummarizedExperiment <- SummarizedExperiment(assays = list(tpm = tpm))

# get s13 phenotype info from column labels
phenoDataS13 <- colnames(tpm)
table(gsub(".*_", "", phenoDataS13))
phenoDataS13 <- data.frame(sample.id = phenoDataS13,
                           source.id = gsub("_.*", "", phenoDataS13),
                           sample.type = gsub(".*_", "", phenoDataS13))
phenoDataS13$tissue.type <- 
  ifelse(phenoDataS13$sample.type=="PBMC", "PBMC", "immune_cell")
phenoDataS13$tissue.type.detail <- 
  ifelse(phenoDataS13$sample.type=="PBMC", "PBMC", phenoDataS13$sample.type)

# append summary statistics to pheno data
tpm <- as.matrix(tpm)
phenoDataS13$library.size <- colSums(tpm)
phenoDataS13$mean.expression <- colMeans(tpm)
phenoDataS13$median.expression <- colMedians(tpm)
phenoDataS13$sd.expression <- colSds(tpm)
phenoDataS13$num.na.expression <- colAnyNAs(tpm)
phenoDataS13$num.zero.expression <- unlist(apply(tpm, 2, function(ci){length(ci[ci==0])}))
rownames(phenoDataS13) <- phenoDataS13[,1]

# append pheno data to colData for SummarizedExperiment
colData(newSummarizedExperiment) <- DataFrame(phenoDataS13)


# map gene symbols using biomaRt
# begin rowdata
newRowData <- data.frame(ensembl_gene_id_version = rownames(tpm))
rownames(newRowData) <- newRowData[,1]
rowData(newSummarizedExperiment) <- DataFrame(newRowData)
# map gene ids to symbols
ensemblMart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
geneIdVector <- rownames(tpm)
geneIdVector <- gsub("\\..*", "", geneIdVector)
rowDataMaps <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = as.character(geneIdVector), 
  mart = ensemblMart
)

# map symbols for rowdata
rowDataSymbol <- rowDataMaps[,2]
names(rowDataSymbol) <- rowDataMaps[,1]
rowDataSymbol <- rowDataSymbol[gsub("\\..*", "", rownames(newSummarizedExperiment))]
length(rowDataSymbol)
dim(newSummarizedExperiment)
rowData(newSummarizedExperiment)$gene_symbol <- rowDataSymbol

#-----
# save
#-----

# save new SummarizedExperiment
save(newSummarizedExperiment, file = "outputs/01_tpm_summaries/se_tpm_s13.rda")

# save env
save.image("./env/01_tpm_summaries/01_read_script.RData")
