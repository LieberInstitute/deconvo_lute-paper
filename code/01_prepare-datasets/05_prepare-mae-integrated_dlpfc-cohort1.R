#!/usr/bin/env R

#
# Makes a multi assay experiment object containing the data for integrated analysis. 
# The new MAE object includes the following asssays:
# * snRNAseq
# * bulk RNAseq
# * RNAscope image processing outputs from HALO
#

libv <- c("here", "nlme", "ggplot2", "gridExtra", "dplyr", "ggforce", "MultiAssayExperiment",
          "SPIAT")
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
# isolate sce data
sce1 <- lscef[[1]]
sce2 <- lscef[[2]]
sce3 <- lscef[[3]]
rm(lscef)

# make sce with image data
img <- halo.output.table
new.img.colnames <- paste0("cell", seq(nrow(img)))
img.data.colnames <- c("Nucleus_Area", "DAPI_Nucleus_Intensity")
img.list <- lapply(img.colnames, function(colname){
  new.data <- img[,colname] %>% as.matrix() %>% t()
  colnames(new.data) <- new.img.colnames
  new.data
})
names(img.list) <- img.colnames
sce.image <- SingleCellExperiment(assays = img.list)
colData(sce.image) <- img

nucleus.area <- img$Nucleus_Area
dapi <- img$DAPI_Nucleus_Intensity


gene.colnames <- c("")
img.expr <- do.call(rbind, lapply())

sce.image <- SingleCellExperiment(assays = list(image = t(halo.output.table)))

# get transpose of large image file
image <- halo.output.table %>% as.data.table() %>% transpose()

colnames(image) <- paste0("cell", seq(ncol(image)))
rm(halo.output.table)

# format image data as very wide data.table
img <- halo.output.table %>% t() %>% as.data.table()
rm(halo.output.table)
colnames(img) <- paste0("cell", seq(ncol(img)))

gc()

#---------
# make mae
#---------
# common sample ids
sample.id.snrnaseq <- "Sample"
sample.id.halo <- "Sample"
sample.id.bulk <- "batch.id2"

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
image.map <- data.frame(colname = colnames(img),
                        primary = img[sample.id.halo,])
listmap <- list(sn.rnaseq = sn.map,
                bulk.rnaseq = bulk.map,
                rnascope.image = image.map)
dfmap <- listToMap(listmap) # make new sampleMap object
rownames(dfmap) <- dfmap$primary

# get object list
object.list <- list(bulk.rnaseq = assays(rse.filter)[["counts"]],
                    sn.rnaseq = assays(sce)[["logcounts"]],
                    rnascope.image = halo.output.table)

# get coldata (harmonized sample ids)
coldata <- data.frame(sample.id = unique(c(sn.map$sample.id, bulk.map$sample.id, image.map$sample.id)))

# make new mae object
mae <- MultiAssayExperiment(experiments = object.list, sampleMap = dfmap, colData = coldata)

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







