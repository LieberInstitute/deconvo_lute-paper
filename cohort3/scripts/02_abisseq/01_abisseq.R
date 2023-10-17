#!/usr/bin/env R

# Author: Sean Maden
#
# Performs deconvolution with the ABIS-seq signature matrix reference.
#
#
#

libv <- c("lute")
sapply(libv, library, character.only = TRUE)

#-----
# load
#-----
zref <- read.csv("./data/zref_abisseq.csv")
ptrue <- read.csv("./data/ptrue_abis.csv")
load("./env/01_tpm_summaries/01_read_script.RData")

#--------------
# format inputs
#--------------
# zref
zref <- zref[!zref[,1]=="12:00 AM",]
rownames(zref) <- zref[,1]
zref <- zref[,c(2:ncol(zref))]

# se
rownames(se) <- rowData(se)$gene_symbol
filter.rows <- !(duplicated(rownames(se))|is.na(rownames(se)))
se <- se[filter.rows,]

#-----------------------
# get cell scale factors
#-----------------------


# filter se marker genes
se <- se[rownames(se) %in% rownames(zref),]
dim(se)

# cell type label mappings
vector.sample.types <- unique(se$sample.type)
dfmap <- data.frame(
  term1 = vector.sample.types,
  term2 = c("B.Naive", "CM", "B.Memory", "")
)

# mean library sizes
s.vector <- unlist(apply(vector.sample.types, function(){
  mean(colSums(
  ))
}))

#------------
# run deconvo
#------------

list.result1 <- lute(
  z = as.matrix(zref), 
  y = as.matrix(assays(se)[["tpm"]]), 
  assay.name = 'tpm',
  typemarker.algorithm = NULL
)

list.result2 <- lute(
  z = as.matrix(zref), 
  y.se = se, 
  celltype.variable = ,
  assay.name = 'tpm',
  typemarker.algorithm = NULL
)



#-----
# save
#-----
save.image("./env/02_abisseq/01_abisseq_script.RData")
