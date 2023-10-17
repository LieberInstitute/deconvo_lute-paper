#!/usr/bin/env R

# Author: Sean Maden
#
# Read in TPM expression data.

libv <- c("SummarizedExperiment", "biomaRt")
sapply(libv, library, character.only = TRUE)

#-----
# load
#-----
tpm.path <- "./data/GSE107011_Processed_data_TPM/GSE107011_Processed_data_TPM.txt"
tpm <- read.table(tpm.path)

#--------------------
# format as se object
#--------------------
se <- SummarizedExperiment(assays = list(tpm = tpm))

#------------------
# get pheno/coldata
#------------------
cd <- colnames(tpm)
table(gsub(".*_", "", cd))
cd <- data.frame(sample.id = cd, 
                 source.id = gsub("_.*", "", cd),
                 sample.type = gsub(".*_", "", cd))
cd$tissue.type <- ifelse(cd$sample.type=="PBMC", "PBMC", "immune_cell")
cd$tissue.type.detail <- ifelse(cd$sample.type=="PBMC", "PBMC", cd$sample.type)
# append summary statistics
tpm <- as.matrix(tpm)
cd$library.size <- colSums(tpm)
cd$mean.expression <- colMeans(tpm)
cd$median.expression <- colMedians(tpm)
cd$sd.expression <- colSds(tpm)
cd$num.na.expression <- colAnyNAs(tpm)
cd$num.zero.expression <- unlist(apply(tpm, 2, function(ci){length(ci[ci==0])}))
rownames(cd) <- cd[,1]

# format summarized experiment
colData(se) <- DataFrame(cd)

#-----------------
# map gene symbols
#-----------------
# begin rowdata
rd.new <- data.frame(ensembl_gene_id_version = rownames(tpm))
rownames(rd.new) <- rd.new[,1]
rowData(se) <- DataFrame(rd.new)

# map gene ids to symbols
mart.ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene.id.vector <- rownames(tpm)
gene.id.vector <- gsub("\\..*", "", gene.id.vector)
rd.maps <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = as.character(gene.id.vector), 
  mart = mart.ensembl
)

# map symbols for rowdata
rd.symbol <- rd.maps[,2]
names(rd.symbol) <- rd.maps[,1]
rd.symbol <- rd.symbol[gsub("\\..*", "", rownames(se))]
length(rd.symbol)
dim(se)
rowData(se)$gene_symbol <- rd.symbol

#-----
# save
#-----
save.image("./env/01_tpm_summaries/01_read_script.RData")