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

# set new mae filename
# new.mae.filename <- "mae_with-rpkm_additional-data_final.rda"
new.mae.filename <- "mae_y-rpkm_z-validate_cohort1.rda"

#--------------------
# load prepped assays
#--------------------
# load snrnaseq data
load("~/GitHub/deconvo_method-paper/outputs/01_prepare-datasets/list-scef_markers-k2-k3-k4_ro1-dlpfc.rda")
load("~/GitHub/deconvo_method-paper/outputs/01_prepare-datasets/list-sce-validate_markers-k-2-3-4_cohort1.rda")
# load bulk data
load("~/GitHub/deconvo_method-paper/outputs/01_prepare-datasets/rse-gene-filter.rda")
load("~/GitHub/deconvo_method-paper/outputs/01_prepare-datasets/rse-rpkmCounts_Human_DLPFC_Deconvolution_n113.rda")
# load rnascope image data
load("~/GitHub/deconvo_method-paper/outputs/01_prepare-datasets/halo-outputs_updated.Rdata")

# common sample ids
sample.id.snrnaseq <- "Sample"
sample.id.halo <- "Sample"
sample.id.bulk <- "batch.id2"

# helper functions
make_sce_all <- function(sce.train, sce.validate){
  # rowdata
  shared.markers <- intersect(rownames(sce.train), rownames(sce.validate))
  sce.train <- sce.train[rownames(sce.train) %in% shared.markers,]
  sce.validate <- sce.validate[rownames(sce.validate) %in% shared.markers,]
  sce.train <- sce.train[!duplicated(rownames(sce.train)),]
  sce.validate <- sce.validate[!duplicated(rownames(sce.validate)),]
  sce.train<- sce.train[order(match(rownames(sce.train), shared.markers)),]
  sce.validate <- sce.validate[order(match(rownames(sce.validate), shared.markers)),]
  # coldata
  cd.sce <- colData(sce.train)
  cd.sce.validate <- colData(sce.validate)
  cn.cd.sce <- colnames(cd.sce)
  cn.cd.sce.validate <- colnames(cd.sce.validate)
  cd.keep <- intersect(cn.cd.sce, cn.cd.sce.validate)
  colData(sce.train) <- cd.sce[,colnames(cd.sce) %in% cd.keep]
  colData(sce.validate) <- cd.sce.validate[,colnames(cd.sce.validate) %in% cd.keep]
  # assays
  assays(sce.train) <- list(counts = assays(sce.train)[["counts"]])
  assays(sce.validate) <- list(counts = assays(sce.validate)[["counts"]])
  rowData(sce.train) <- rowData(sce.validate)
  metadata(sce.train) <- metadata(sce.validate) <- list()
  # checks
  identical(names(assays(sce.train)), names(assays(sce.validate)))
  identical(rowData(sce.train), rowData(sce.validate))
  identical(rownames(sce.train), rownames(sce.validate))
  identical(colnames(colData(sce.train)), colnames(colData(sce.validate)))
  # make sce.all
  metadata <- list(train = metadata(sce.train),
                   validate = metadata(sce.validate))
  assays.all <- list(counts = cbind(counts(sce.train), counts(sce.validate)))
  sce.all <- SingleCellExperiment(assays = assays.all,
                                  rowData = rowData(sce.train))
  colData(sce.all) <- rbind(colData(sce.train), colData(sce.validate))
  metadata(sce.all) <- metadata
  return(sce.all)
}

#------------------------------------
# snrnaseq: isolate snrnaseq sce data
#------------------------------------
# train
sce1 <- lscef[[1]]
sce2 <- lscef[[2]]
sce3 <- lscef[[3]]
rm(lscef)
# validate
lscef.validate <- list.sce.validate.markers
sce1.validate <- lscef.validate[[1]]
sce2.validate <- lscef.validate[[2]]
sce3.validate <- lscef.validate[[3]]
rm(lscef.validate)
rm(list.sce.validate.markers)

# all (train + validate)
sce1.all <- make_sce_all(sce1, sce1.validate)
sce2.all <- make_sce_all(sce2, sce2.validate)
sce3.all <- make_sce_all(sce3, sce3.validate)

#--------------------------------------------
# bulk: isolate bulk rnaseq marker expression
#--------------------------------------------
# match counts and rpkm data
rownames(rse.filter) <- rowData(rse.filter)$Symbol
rownames(rse.rpkm) <- rowData(rse.rpkm)$Symbol
rse.rpkm <- rse.rpkm[rownames(rse.rpkm) %in% rownames(rse.filter),]
dim(rse.rpkm)

# match coldata
colnames.counts <- colnames(rse.filter)
colnames.counts.filter <- !grepl("_Cyto|_Nuc", colnames(rse.filter))
colnames.counts[colnames.counts.filter] <- paste0(colnames.counts[colnames.counts.filter],"_Bulk")
colnames(rse.filter) <- colnames.counts
head(colnames(rse.filter))
head(rownames(colData(rse.filter)))
identical(colnames(rse.filter), colnames(rse.rpkm))
colData(rse.rpkm) <- colData(rse.filter)

# filter markers
marker.id.vector <- unique(c(rownames(sce1), rownames(sce2), rownames(sce3)))
rse.filter <- rse.filter[rownames(rse.filter) %in% marker.id.vector,]
rse.rpkm <- rse.rpkm[rownames(rse.rpkm) %in% marker.id.vector,]

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
# make transpose
df.rnascope.kdata <- t(df.rnascope.kdata)

#-------------
# mae: coldata
#-------------
# coldata
unique.sample.id.vector <- unique(
  c(sce1[[sample.id.snrnaseq]], 
    sce1.validate[[sample.id.snrnaseq]],
    rse.filter[[sample.id.bulk]], # bulk and rpkm have same sample ids ...
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
sn1.validate.map <- data.frame(colname = colnames(sce1.validate), 
                               primary = sce1.validate[[sample.id.snrnaseq]])
sn2.validate.map <- data.frame(colname = colnames(sce2.validate), 
                               primary = sce2.validate[[sample.id.snrnaseq]])
sn3.validate.map <- data.frame(colname = colnames(sce3.validate), 
                               primary = sce3.validate[[sample.id.snrnaseq]])
sn1.all.map <- data.frame(colname = colnames(sce1.all), 
                               primary = sce1.all[[sample.id.snrnaseq]])
sn2.all.map <- data.frame(colname = colnames(sce2.all), 
                               primary = sce2.all[[sample.id.snrnaseq]])
sn3.all.map <- data.frame(colname = colnames(sce3.all), 
                               primary = sce3.all[[sample.id.snrnaseq]])
# bulk
bulk.map <- data.frame(colname = colnames(rse.filter), primary = rse.filter[[sample.id.bulk]])
bulk.rpkm.map <- data.frame(colname = colnames(rse.rpkm), primary = rse.rpkm[[sample.id.bulk]])
# rnascope data
image.map <- data.frame(colname = colnames(sce.img), primary = sce.img[[sample.id.halo]])
dfrn.map <- data.frame(colname = colnames(df.rnascope.kdata), 
                       primary = gsub(";.*", "", colnames(df.rnascope.kdata)))

# make mappings list, then map df
listmap <- list(sn1.rnaseq = sn1.map, # snrnaseq
                sn2.rnaseq = sn2.map, 
                sn3.rnaseq = sn3.map,
                sn1.validate.rnaseq = sn1.validate.map, # snrnaseq
                sn2.validate.rnaseq = sn2.validate.map, 
                sn3.validate.rnaseq = sn3.validate.map,
                sn1.all.rnaseq = sn1.all.map, # snrnaseq
                sn2.all.rnaseq = sn2.all.map, 
                sn3.all.rnaseq = sn3.all.map,
                bulk.rnaseq = bulk.map, # bulk
                bulk.rpkm.rnaseq = bulk.rpkm.map, # bulk, rpkm
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
  sn1.validate.rnaseq = sce1.validate, 
  sn2.validate.rnaseq = sce2.validate, 
  sn3.validate.rnaseq = sce3.validate, 
  sn1.all.rnaseq = sce1.all, 
  sn2.all.rnaseq = sce2.all, 
  sn3.all.rnaseq = sce3.all, 
  bulk.rnaseq = rse.filter,
  bulk.rpkm.rnaseq = rse.rpkm,
  rnascope.image = sce.img,
  df.cellstat.rnascope = df.rnascope.kdata %>% as.matrix())
experiment.list <- ExperimentList(object.list)
rownames(coldata) <- coldata$sample.id

#------------------------------
# mae: prepare new mae datasets
#------------------------------
mae <- prepMultiAssay(ExperimentList = experiment.list, 
                      sampleMap = dfmap, colData = coldata)
mae.final <- MultiAssayExperiment(mae$experiments, 
                                  mae$colData, 
                                  mae$sampleMap)

#----------
# mae: save
#----------
mae.final.filepath <- here("outputs", "01_prepare-datasets", new.mae.filename)
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