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
# load snrnaseq validation data
load("~/GitHub/deconvo_method-paper/outputs/01_prepare-datasets/list-sce-validate_k-2-3-4-markers_cohort1.rda")
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

# isolate snrnaseq sce data
sce1 <- lscef[[1]]
sce2 <- lscef[[2]]
sce3 <- lscef[[3]]
rm(lscef)
# isolate snrnaseq sce data -- VALIDATION
lscef.validate <- list.sce.validate.markers
sce1.validate <- lscef.validate[[1]]
sce2.validate <- lscef.validate[[2]]
sce3.validate <- lscef.validate[[3]]
rm(list.sce.validate.markers)
rm(lscef.validate)

# make sce with image data
img <- halo.output.table
new.img.colnames <- paste0("cell", seq(nrow(img)))
img.data.colnames <- c("Nucleus_Area", "AKT3_Copies", "Cell_Area", 
                       "DAPI_Nucleus_Intensity", "DAPI_Cytoplasm_Intensity")
img.coldata.colnames <- c("SAMPLE_ID", sample.id.halo, "cell_type", "Slide", 
                          "XMin", "XMax", "YMin", "YMax")
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
# indicies for subsetting
bulk.index <- seq(ncol(rse.filter))
sn.index <- seq(ncol(sce1))
sn.validate.index <- seq(ncol(sce1.validate))
img.index <- seq(ncol(sce.img))

# mae prep

# get subsets
bulk.rse.subset <- rse.filter[,bulk.index]
sn.sce.subset <- sce1[,sn.index]
sn.sce.validate.subset <- sce1.validate[,sn.validate.index]
img.sce.subset <- sce.img[,img.index]

# coldata
coldata <- DataFrame(data.frame(sample.id = unique(c(sn.sce.subset[[sample.id.snrnaseq]],
                                                     sn.sce.validate.subset[[sample.id.snrnaseq]],
                                                     bulk.rse.subset[[sample.id.bulk]],
                                                     img.sce.subset[[sample.id.halo]]))))

# make maplist
sn1.map <- data.frame(colname = colnames(sn.sce.subset), 
                      primary = sn.sce.subset[[sample.id.snrnaseq]])
sn1.validate.map <- data.frame(colname = colnames(sn.sce.validate.subset), 
                      primary = sn.sce.validate.subset[[sample.id.snrnaseq]])
bulk.map <- data.frame(colname = colnames(bulk.rse.subset), 
                       primary = bulk.rse.subset[[sample.id.bulk]])
image.map <- data.frame(colname = colnames(img.sce.subset), 
                        primary = img.sce.subset[[sample.id.halo]])
listmap <- list(sn1.rnaseq = sn1.map,
                sn1.validate.rnaseq = sn1.validate.map,
                bulk.rnaseq = bulk.map, 
                rnascope.image = image.map)
dfmap <- listToMap(listmap) # make new sampleMap object
rownames(dfmap) <- dfmap$primary
object.list2 <- list(bulk.rnaseq = bulk.rse.subset, 
                     sn1.rnaseq = sn.sce.subset, 
                     sn1.validate.rnaseq = sn.sce.validate.subset, 
                     rnascope.image = img.sce.subset)
experiment.list <- ExperimentList(object.list2)
rownames(coldata) <- coldata$sample.id
mae <- prepMultiAssay(ExperimentList = experiment.list, sampleMap = dfmap, colData = coldata)
mae.final <- MultiAssayExperiment(mae$experiments, mae$colData, mae$sampleMap)

# save
mae.final.filepath <- here("deconvo_method-paper", "outputs", "01_prepare-datasets", "mae_final.rda")
save(mae.final, file = mae.final.filepath)

#--------------
# mae summaries
#--------------
# summarize
complete.id <- colData(mae.final)[complete.cases(mae.final),]
length(complete.id) # 11

# upset plot of samples by assays
upsetSamples(mae.final)

# example: subset on one sample id
mae.final[,colData(mae.final)$sample.id=="Br8492_mid",]

#-----------------------------------
# mae: inspect, with basic summaries
#-----------------------------------
experiments(mae)
colData(mae)
