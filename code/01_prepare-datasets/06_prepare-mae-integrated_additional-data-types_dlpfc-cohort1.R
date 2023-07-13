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

# common sample ids
sample.id.snrnaseq <- "Sample"
sample.id.halo <- "Sample"
sample.id.bulk <- "batch.id2"

#------------------------------------
# snrnaseq: isolate snrnaseq sce data
#------------------------------------
sce1 <- lscef[[1]]
sce2 <- lscef[[2]]
sce3 <- lscef[[3]]
rm(lscef)

#--------------------------------------------
# bulk: isolate bulk rnaseq marker expression
#--------------------------------------------
marker.id.vector <- unique(c(rownames(sce1), rownames(sce2), rownames(sce3)))
rownames(rse.filter) <- rowData(rse.filter)$Symbol
rse.filter <- rse.filter[rownames(rse.filter) %in% marker.id.vector,]

#-----------------------------------
# rnascope: make sce with image data
#-----------------------------------
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

#--------------------------------------------------
# df.rn: get additional rnascope data.frame objects
#--------------------------------------------------
# update sce.img with klabel categories
# get true proportions from rnascope data
rnascope <- sce.img
cd <- colData(rnascope)
cd$k2 <- cd$k3 <- cd$k4 <- "NA"
# append k label categories
cd$k2 <- ifelse(grepl("Excit|Inhib", cd$cell_type), "neuron",
                ifelse(grepl("Endo|Oligo|Micro", cd$cell_type), "glial", "NA"))
cd$k3 <- ifelse(grepl("Excit", cd$cell_type), "Excit",
                ifelse(grepl("Inhib", cd$cell_type), "Inhib", 
                       ifelse(grepl("Endo|Oligo|Micro", cd$cell_type), "glial", "NA")))
cd$k4 <- ifelse(grepl("Excit", cd$cell_type), "Excit",
                ifelse(grepl("Inhib", cd$cell_type), "Inhib", 
                       ifelse(grepl("Oligo", cd$cell_type), "Oligo", 
                              ifelse(grepl("Micro|Endo", cd$cell_type), "non_oligo_glial", "NA"))))
# reassign coldata to sce.img
colData(sce.img) <- cd

# get data.frames of cell sizes, counts, proportions by klabel category
sample.id.vector <- unique(rnascope$Sample)
# filter na values
rnascope <- sce.img[,!sce.img$k2=="NA"]
# get all kdata
df.rnascope.kdata <- do.call(rbind, lapply(c("k2", "k3", "k4"), function(cell.type.variable){
  do.call(rbind, lapply(sample.id.vector, function(sample.id){
    rnf <- rnascope[,rnascope$Sample==sample.id]
    # proportions
    df.prop <- table(rnf[[cell.type.variable]], rnf$Sample) %>% prop.table()
    # counts
    df.count <- table(rnf[[cell.type.variable]], rnf$Sample)
    # sizes
    df.size <- aggregate(data.frame(area = assays(rnf)[["Nucleus_Area"]][1,]), 
                         by = list(cell_type = rnf[[cell.type.variable]]), FUN = "median")
    # combine and format output
    df.iter <- cbind(df.prop, cbind(df.size, df.count))
    df.iter <- df.iter[,c(1,2,3,5,8)]
    colnames(df.iter) <- c("cell_type", "sample_id", "true_proportion", "cell_size", "cell_count")
    df.iter$k.label <- cell.type.variable
    df.iter
  }))
}))
rownames(df.rnascope.kdata) <- paste0(df.rnascope.kdata$sample_id,";",
                                      df.rnascope.kdata$cell_type, ";",
                                      df.rnascope.kdata$k.label)


#-------------
# mae: coldata
#-------------
# coldata
unique.sample.id.vector <- unique(
  c(sce1[[sample.id.snrnaseq]], 
    rse.filter[[sample.id.bulk]],
    sce.img[[sample.id.halo]])
)
coldata <- DataFrame(
  data.frame(
    sample.id = unique.sample.id.vector))

#------------------------
# mae: make data mappings
#------------------------
# make maplist
# snrnaseq
sn1.map <- data.frame(colname = colnames(sce1), primary = sce1[[sample.id.snrnaseq]])
sn2.map <- data.frame(colname = colnames(sce2), primary = sce2[[sample.id.snrnaseq]])
sn3.map <- data.frame(colname = colnames(sce3), primary = sce3[[sample.id.snrnaseq]])
# bulk
bulk.map <- data.frame(colname = colnames(rse.filter), primary = rse.filter[[sample.id.bulk]])
# rnascope data
image.map <- data.frame(colname = colnames(sce.img), primary = sce.img[[sample.id.halo]])
dfrn.map <- data.frame(colname = rownames(df.rnascope.kdata), primary = df.rnascope.kdata$sample_id)

# make mappings list, then map df
listmap <- list(sn1.rnaseq = sn1.map, # snrnaseq
                sn2.rnaseq = sn2.map, 
                sn3.rnaseq = sn3.map,
                bulk.rnaseq = bulk.map, # bulk
                rnascope.image = image.map, # rnascope
                df.cellstat.rnascope = dfrn.map)
dfmap <- listToMap(listmap) # make new sampleMap object
rownames(dfmap) <- dfmap$primary

#---------------------
# mae: get object list
#---------------------
object.list <- list(
  sn1.rnaseq = sce1, 
  sn2.rnaseq = sce2, 
  sn3.rnaseq = sce3, 
  bulk.rnaseq = rse.filter,
  rnascope.image = sce.img,
  df.cellstat.rnascope = df.rnascope.kdata)   
experiment.list <- ExperimentList(object.list)
rownames(coldata) <- coldata$sample.id

#------------------------------
# mae: prepare new mae datasets
#------------------------------
mae <- prepMultiAssay(ExperimentList = experiment.list, sampleMap = dfmap, colData = coldata)
mae.final <- MultiAssayExperiment(mae$experiments, mae$colData, mae$sampleMap)

#----------
# mae: save
#----------
mae.final.filepath <- here("deconvo_method-paper", "outputs", "01_prepare-datasets", "mae_additional-data_final.rda")
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
