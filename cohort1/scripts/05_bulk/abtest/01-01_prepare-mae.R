#!/usr/bin/env R

# Prepares MultiAssayExperiment object

# experiment variables
deconvolution.algorithm <- "nnls"
cell.type.variable <- "k2"

#-----------
# unpack mae
#-----------
# get bulk expression
rse.counts <- mae[["bulk.rnaseq"]]
rse.rpkm <- mae[["bulk.rpkm.rnaseq"]]
names(assays(rse.counts)) <- names(assays(rse.rpkm)) <- "counts"

# snrnaseq reference -- using same reference across experiments
sce.iter <- mae[["snrnaseq.k2.all"]]
sce.iter <- logNormCounts(sce.iter)

# prep rnascope data
# get true proportions from rnascope data
df.rn <- mae[["cell.sizes"]]
df.rn <- as.data.frame(t(df.rn))
df.rn <- df.rn[df.rn$k.label=="k2",]

# get list.df.true
list.df.true <- metadata(sce.iter)[["list.df.true.k2"]]
