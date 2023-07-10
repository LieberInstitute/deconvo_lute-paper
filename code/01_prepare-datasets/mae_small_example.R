
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

---------
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
object.list <- list(bulk.rnaseq = assays(rse.filter)[["counts"]],
                    sn.rnaseq = assays(sce)[["logcounts"]],
                    rnascope.image = halo.output.table)

experiment.list <- ExperimentList(object.list)