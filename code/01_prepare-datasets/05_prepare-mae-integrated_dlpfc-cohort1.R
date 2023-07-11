#!/usr/bin/env R

#
# Makes a multi assay experiment object containing the data for integrated analysis. 
# The new MAE object includes the following asssays:
# * snRNAseq
# * bulk RNAseq
# * RNAscope image processing outputs from HALO
#

libv <- c("here", "nlme", "ggplot2", "gridExtra", "dplyr", "ggforce", "MultiAssayExperiment", "SingleCellExperiment")
sapply(libv, library, character.only = TRUE)

#--------------------
# load prepped assays
#--------------------
# load snrnaseq data
load("~/GitHub/deconvo_method-paper/outputs/01_prepare-datasets/list-scef_markers-k2-k3-k4_ro1-dlpfc.rda")
# load bulk data
load("~/GitHub/deconvo_method-paper/outputs/01_prepare-datasets/rse-gene-filter.rda")
# load rnascope image data
load("~/GitHub/deconvo_method-paper/outputs/01_prepare-datasets/halo-outputs_updated.Rdata")

#------------------
# prep data for mae
#------------------
# common sample ids
sample.id.snrnaseq <- "Sample"
sample.id.halo <- "Sample"
sample.id.bulk <- "batch.id2"

# isolate sce data
sce1 <- lscef[[1]]
sce2 <- lscef[[2]]
sce3 <- lscef[[3]]
rm(lscef)

# make sce with image data
img <- halo.output.table
new.img.colnames <- paste0("cell", seq(nrow(img)))
img.data.colnames <- c("Nucleus_Area", "AKT3_Copies", "Cell_Area", 
                       "DAPI_Nucleus_Intensity", "DAPI_Cytoplasm_Intensity")
img.coldata.colnames <- c("SAMPLE_ID", sample.id.halo, "Slide", "XMin", "XMax", "YMin", "YMax")
img.list <- lapply(img.data.colnames, function(colname){
  new.data <- img[,colname] %>% as.matrix() %>% t()
  colnames(new.data) <- new.img.colnames
  new.data
})
names(img.list) <- img.data.colnames
sce.img <- SingleCellExperiment(assays = img.list)
img.coldata <- DataFrame(img[,img.coldata.colnames])
rownames(img.coldata) <- new.img.colnames
colData(sce.img) <- img.coldata
rm(halo.output.table)

gc()

#---------
# make mae
#---------
# make mae on common indices
# make maplist
sn1.map <- data.frame(colname = colnames(sce1),
                     primary = sce1[[sample.id.snrnaseq]])
sn2.map <- data.frame(colname = colnames(sce2),
                      primary = sce2[[sample.id.snrnaseq]])
sn3.map <- data.frame(colname = colnames(sce3),
                      primary = sce3[[sample.id.snrnaseq]])
bulk.map <- data.frame(colname = colnames(rse.filter),
                       primary = rse.filter[[sample.id.bulk]])
image.map <- data.frame(colname = colnames(sce.img),
                        primary = sce.img[[sample.id.halo]])
listmap <- list(sn1.rnaseq = sn1.map,
                sn2.rnaseq = sn2.map,
                sn3.rnaseq = sn3.map,
                bulk.rnaseq = bulk.map,
                rnascope.image = image.map)
dfmap <- listToMap(listmap) # make new sampleMap object
rownames(dfmap) <- dfmap$primary

# get object list
object.list <- list(bulk.rnaseq = assays(rse.filter)[["counts"]],
                    sn1.rnaseq = assays(sce1)[["logcounts"]],
                    sn2.rnaseq = assays(sce2)[["logcounts"]],
                    sn3.rnaseq = assays(sce3)[["logcounts"]],
                    rnascope.image = sce.img)

# get coldata (harmonized sample ids)
coldata <- data.frame(sample.id = unique(c(sn1.map$sample.id,
                                           sn2.map$sample.id,
                                           sn3.map$sample.id, 
                                           bulk.map$sample.id, 
                                           image.map$sample.id)))

# make new mae object
mae <- MultiAssayExperiment(experiments = object.list, 
                            sampleMap = dfmap, 
                            colData = coldata)

###


experiment.list <- ExperimentList(object.list)

###

mae <- prepMultiAssay(ExperimentList = experiment.list, sampleMap = dfmap, colData = coldata)
experiments(mae)

mae <- prepMultiAssay(ExperimentList = experiment.list)
experiments(mae)

mae <- prepMultiAssay(ExperimentList = experiment.list, colData = coldata)
experiments(mae)

mae <- prepMultiAssay(ExperimentList = experiment.list, sampleMap = dfmap)
experiments(mae)

#-----------------------------------
# mae: inspect, with basic summaries
#-----------------------------------
experiments(mae)

colData(mae)

mae <- prepMultiAssay(ExperimentList = ExperimentList(bulk.rnaseq = rse.filter), 
                      sampleMap = dfmap[dfmap$assay=="bulk.rnaseq",])


