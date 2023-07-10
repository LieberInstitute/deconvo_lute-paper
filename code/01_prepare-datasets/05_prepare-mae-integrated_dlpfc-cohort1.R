#!/usr/bin/env R

#
# Makes a multi assay experiment object containing the data for integrated analysis. 
# The new MAE object includes the following asssays:
# * snRNAseq
# * bulk RNAseq
# * RNAscope image processing outputs from HALO
#

libv <- c("here", "nlme", "ggplot2", "gridExtra", "dplyr", "ggforce", "MultiAssayExperiment")
sapply(libv, library, character.only = TRUE)

#--------------------
# load prepped assays
#--------------------
# load snrnaseq data
load("~/GitHub/deconvo_method-paper/outputs/01_prepare-datasets/list-scef_markers-k2-k3-k4_ro1-dlpfc.rda")
sce <- lscef[[1]]
# load bulk data
load("~/GitHub/deconvo_method-paper/outputs/01_prepare-datasets/rse-gene-filter.rda")
# load rnascope image data
load("~/GitHub/deconvo_method-paper/outputs/01_prepare-datasets/halo-outputs_updated.Rdata")

#---------
# make mae
#---------
sample.id.snrnaseq <- "Sample"
sample.id.halo <- "Sample"
sample.id.bulk <- "batch.id2"

# make mae on common indices
# make maplist
sn.map <- data.frame(coldata = colnames(sce),
                     sample.id = sce[[sample.id.snrnaseq]])
bulk.map <- data.frame(coldata = colnames(rse.filter),
                       sample.id = rse.filter[[sample.id.bulk]])
image.map <- data.frame(coldata = halo.output.table[,sample.id.halo],
                        sample.id = halo.output.table[,sample.id.halo])
listmap <- list(sn.rnaseq = sn.map,
                bulk.rnaseq = bulk.map,
                rnascope.image = image.map)

# make new sampleMap object
listToMap(maplist)


