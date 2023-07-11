
libv <- c("here", "nlme", "ggplot2", "gridExtra", "dplyr", "ggforce", "MultiAssayExperiment")
sapply(libv, library, character.only = TRUE)

#--------------------
# load prepped assays
#--------------------
# load snrnaseq data
load("~/GitHub/deconvo_method-paper/outputs/01_prepare-datasets/list-scef_markers-k2-k3-k4_ro1-dlpfc.rda")
sce <- lscef[[1]][seq(10),]
rm(lscef)

# load bulk data
load("~/GitHub/deconvo_method-paper/outputs/01_prepare-datasets/rse-gene-filter.rda")
rse <- rse.filter[seq(10),]
rm(rse.filter)

# load rnascope image data
load("~/GitHub/deconvo_method-paper/outputs/01_prepare-datasets/halo-outputs_updated.Rdata")
image <- halo.output.table[seq(10),]
rownames(image) <- paste0("cell", seq(nrow(image)))
image <- t(image)
rm(halo.output.table)
gc()

sample.id.snrnaseq <- "Sample"
sample.id.halo <- "Sample"
sample.id.bulk <- "batch.id2"

# make mae on common indices
# make maplist
sn.map <- data.frame(colname = colnames(sce),
                     primary = sce[[sample.id.snrnaseq]])
bulk.map <- data.frame(colname = colnames(rse),
                       primary = rse[[sample.id.bulk]])
image.map <- data.frame(colname = colnames(image),
                        primary = as.character(image[sample.id.halo,]))
listmap <- list(sn.rnaseq = sn.map,
                bulk.rnaseq = bulk.map,
                rnascope.image = image.map)
dfmap <- listToMap(listmap) # make new sampleMap object
rownames(dfmap) <- dfmap$primary

# get object list
object.list <- list(bulk.rnaseq = assays(rse)[["counts"]],
                    sn.rnaseq = assays(sce)[["counts"]],
                    rnascope.image = as.matrix(image))
experiment.list <- ExperimentList(object.list)
mae <- MultiAssayExperiment(experiments = object.list, 
                            sampleMap = dfmap)
experiments(mae)







# get object list
object.list <- list(bulk.rnaseq = assays(rse)[["counts"]],
                    sn.rnaseq = assays(sce)[["counts"]])
experiment.list <- ExperimentList(object.list)
mae <- MultiAssayExperiment(experiments = object.list, 
                            sampleMap = dfmap)
experiments(mae)

# get object list
object.list <- list(bulk.rnaseq = assays(rse)[["counts"]],
                    rnascope.image = as.matrix(image))
experiment.list <- ExperimentList(object.list)
mae <- MultiAssayExperiment(experiments = object.list, 
                            sampleMap = dfmap)
experiments(mae)

# get object list
object.list <- list(sn.rnaseq = assays(sce)[["counts"]],
                    rnascope.image = as.matrix(image))
experiment.list <- ExperimentList(object.list)
mae <- MultiAssayExperiment(experiments = object.list, 
                            sampleMap = dfmap)
experiments(mae)


