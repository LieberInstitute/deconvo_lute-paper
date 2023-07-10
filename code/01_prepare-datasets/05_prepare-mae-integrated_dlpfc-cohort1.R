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
sn.map <- data.frame(colname = colnames(sce),
                     primary = sce[[sample.id.snrnaseq]])
bulk.map <- data.frame(colname = colnames(rse.filter),
                       primary = rse.filter[[sample.id.bulk]])
image.map <- data.frame(colname = halo.output.table[,sample.id.halo],
                        primary = halo.output.table[,sample.id.halo])
listmap <- list(sn.rnaseq = sn.map,
                bulk.rnaseq = bulk.map,
                rnascope.image = image.map)
dfmap <- listToMap(listmap) # make new sampleMap object
rownames(dfmap) <- dfmap$primary

# get object list
object.list <- list(bulk.rnaseq = rse.filter,
                    sn.rnaseq = sce,
                    rnascope.image = halo.output.table)

# get coldata (harmonized sample ids)
coldata <- data.frame(sample.id = unique(c(sn.map$sample.id, bulk.map$sample.id, image.map$sample.id)))

# make new mae object
mae <- MultiAssayExperiment(experiments = object.list, sampleMap = dfmap, colData = coldata)

#-----------------------------------
# mae: inspect, with basic summaries
#-----------------------------------
experiments(mae)

colData(mae)




