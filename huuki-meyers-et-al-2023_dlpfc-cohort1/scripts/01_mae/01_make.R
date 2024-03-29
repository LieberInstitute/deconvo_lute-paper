#!/usr/bin/env R

# Author: Sean Maden
#
# Run this script from: ./deconvo_method-paper/cohort1/
#
# Makes an object of type MultiAssayExperiment containing experiment data organizead on sample ID 
# (BrNum + Position). Data types in the MultiAssayExperiment:
#
# * snRNAseq
#
# * bulk RNAseq
#
# * RNAscope image processing outputs from HALO
#
#

#-------------
# dependencies
#-------------

libv <- c("here", "nlme", "lute", "ggplot2", "gridExtra", "dplyr", "ggforce", "MultiAssayExperiment", "SingleCellExperiment")
sapply(libv, library, character.only = TRUE)

#--------------------
# load prepped assays
#--------------------
# load snrnaseq data
load("./outputs/00_preprocess/list-scef_markers-k2-k3-k4_ro1-dlpfc.rda")
load("./outputs/00_preprocess/list-sce-validate_k-2-3-4-markers_cohort1.rda")
# load bulk data
load("./outputs/00_preprocess/rse-gene-filter.rda")
load("./outputs/00_preprocess/rse-rpkmCounts_Human_DLPFC_Deconvolution_n113.rda")
# load rnascope image data
load("./data/01_mae/halo_all.Rdata")

#------------------
# common sample ids
#------------------
sample.id.snrnaseq <- "Sample"
sample.id.halo <- "Sample"
sample.id.bulk <- "batch.id2"
sample.id.sn.proportions <- "Sample"

#-----------------
# helper functions
#-----------------
make_sce_all <- function(sce.train, sce.validate){
  # make_sce_all
  #
  # get combined set of train and validation snRNAseq data
  #
  #
  
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

# append cell type category names
celltypevar <- "cellType_broad_hc"
# define marker categories
sce1.validate[["k2"]] <- ifelse(grepl("^Excit.*|^Inhib.*", sce1.validate[[celltypevar]]), 
                      "neuron", "glial")
sce2.validate[["k3"]] <- ifelse(grepl("^Excit.*", sce2.validate[[celltypevar]]), "Excit", 
                      ifelse(grepl("^Inhib.*", sce2.validate[[celltypevar]]), 
                             "Inhib", "glial"))
sce3.validate[["k4"]] <- ifelse(grepl("^Excit.*", sce3.validate[[celltypevar]]), "Excit", 
                      ifelse(grepl("^Inhib.*", sce3.validate[[celltypevar]]), "Inhib", 
                             ifelse(grepl("^Oligo$", sce3.validate[[celltypevar]]), "Oligo", 
                                    "non_oligo_glial")))

# all (train + validate)
sce1.all <- make_sce_all(sce1, sce1.validate)
sce2.all <- make_sce_all(sce2, sce2.validate)
sce3.all <- make_sce_all(sce3, sce3.validate)

#-------------------------------------------
# write the sample ids for train, validation
#-------------------------------------------
sample.id.train <- unique(sce1[[sample.id.snrnaseq]])
sample.id.validate <- unique(sce1.validate[[sample.id.snrnaseq]])
list.sample.id.snrnaseq <- list(train = sample.id.train, validation = sample.id.validate)

# save
save(list.sample.id.snrnaseq, file = "outputs/00_preprocess/list_snrnaseq_sampleid.rda")

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

#-----------
# pseudobulk
#-----------
# make pseudobulk from sce
# for k2, k3, k4
list.sce.k234 <- list(k2 = sce1.all, k3 = sce2.all, k4 = sce3.all)
list.pb.k234 <- lapply(seq(3), function(index){
  sce <- list.sce.k234[[index]]
  sample.id.vector <- unique(sce[["Sample"]])
  variable.name <- paste0("k", index+1)
  # get pseudobulk
  ypb <- do.call(cbind, lapply(sample.id.vector, function(sample.id){
    scef <- sce[,sce[["Sample"]]==sample.id]
    ypb_from_sce(scef, "counts", variable.name)
  }))
  ypb <- as.matrix(ypb)
  colnames(ypb) <- sample.id.vector
  bulk.pb <- SummarizedExperiment(assays = SimpleList(list(counts = ypb)))
  colnames(bulk.pb) <- sample.id.vector
  cd <- DataFrame(data.frame(sample.id = colnames(ypb), batch.id2 = colnames(ypb)))
  rownames(cd) <- colnames(ypb)
  colData(bulk.pb) <- cd
  colnames(bulk.pb)
  as(bulk.pb, "RangedSummarizedExperiment")
})
names(list.pb.k234) <- names(list.sce.k234)

#-----------------------------------
# rnascope: make sce with image data
#-----------------------------------
img <- as.data.frame(as.matrix(halo_all))
new.img.colnames <- paste0("cell", seq(nrow(img)))
img.data.colnames <- c("Nucleus_Area", "AKT3_Copies", "Cell_Area", 
                       "DAPI_Nucleus_Intensity", "DAPI_Cytoplasm_Intensity")
img.coldata.colnames <- c("SAMPLE_ID", sample.id.halo, "cell_type", "Slide", 
                          "XMin", "XMax", "YMin", "YMax", "Confidence")
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
rm(halo_all)
gc()

expected.samples.count.rnascope <- 19
length(unique(sce.img$Sample)) == expected.samples.count.rnascope

#--------------------
# preprocess rnascope
#--------------------
max.nucleus.area <- 78
expected.nuclei.count <- 1362399

# filter on nucleus area
filter.sce <- assays(sce.img)[["Nucleus_Area"]] < max.nucleus.area
sce.img <- sce.img[,filter.sce]

# test: total rnascope nuclei
ncol(sce.img)==expected.nuclei.count

# test: samples overlapping rnascope and snrnaseq
expected.samples.overlapping <- 11
length(intersect(unique(sce.img$Sample), unique(sce1$Sample)))==
  expected.samples.overlapping









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
    df.size <- aggregate(data.frame(area = as.numeric(assays(rnf)[["Nucleus_Area"]][1,])), 
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

# test: number of samples in 1. cell.sizes and 2. sce.img
length(unique(df.rnascope.kdata["sample_id",]))==expected.samples.count.rnascope
identical(unique(df.rnascope.kdata["sample_id",]), unique(sce.img$Sample))

# save rnascope data objects
save(sce.img, file = "outputs/00_preprocess/sce_img.rda")
save(df.rnascope.kdata, file = "outputs/00_preprocess/df_rnascope_kdata.rda")

#-------------
# mae: coldata
#-------------
# coldata
unique.sample.id.vector <- unique(
  c(
    
    sce1[[sample.id.snrnaseq]], 
    
    sce1.validate[[sample.id.snrnaseq]],
    
    rse.filter[[sample.id.bulk]], # bulk and rpkm have same sample ids ...
    
    df.rnascope.kdata["sample_id",],
    
    list.sce.k234[["k2"]][["sample.id"]]
    
    )
)
coldata <- DataFrame(
  data.frame(
    sample.id = unique.sample.id.vector))




#------------------------
# mae: make data mappings/maplist
#------------------------


# snrnaseq

sn.k2.map <- data.frame(colname = colnames(sce1.all), primary = sce1.all[[sample.id.snrnaseq]])

sn.k3.map <- data.frame(colname = colnames(sce2.all), primary = sce2.all[[sample.id.snrnaseq]])

sn.k4.map <- data.frame(colname = colnames(sce3.all), primary = sce3.all[[sample.id.snrnaseq]])


# bulk

bulk.map <- data.frame(colname = colnames(rse.filter), 
                       primary = rse.filter[[sample.id.bulk]])

bulk.rpkm.map <- data.frame(colname = colnames(rse.rpkm), 
                            primary = rse.rpkm[[sample.id.bulk]])

# rnascope data

sce.img.map <- data.frame(colname = colnames(sce.img), 
                        primary = sce.img[[sample.id.halo]])

dfrn.map <- data.frame(colname = colnames(df.rnascope.kdata), 
                       primary = gsub(";.*", "", colnames(df.rnascope.kdata)))

# pseudobulk

bulk.pb.k2.map <- data.frame(colname = colnames(list.pb.k234[["k2"]]),
                             primary = list.pb.k234[["k2"]][["sample.id"]])
bulk.pb.k3.map <- data.frame(colname = colnames(list.pb.k234[["k3"]]),
                             primary = list.pb.k234[["k3"]][["sample.id"]])
bulk.pb.k4.map <- data.frame(colname = colnames(list.pb.k234[["k4"]]),
                             primary = list.pb.k234[["k4"]][["sample.id"]])

# make mappings list, then map df

listmap <- list(
  snrnaseq.k2.all = sn.k2.map, # snrnaseq
                
  snrnaseq.k3.all = sn.k3.map, 
                
  snrnaseq.k4.all = sn.k4.map,
                
  bulk.rnaseq = bulk.map, # bulk
                
  bulk.rpkm.rnaseq = bulk.rpkm.map, # bulk, rpkm
               
  cell.sizes = dfrn.map, # rnascope cell sizes
  
  sce.img = sce.img.map, # rnascope sce image data
  
  bulk.pb.k2 = bulk.pb.k2.map, # pseudobulk k2
  
  bulk.pb.k3 = bulk.pb.k3.map, # pseudobulk k3
  
  bulk.pb.k4 = bulk.pb.k4.map # pseudobulk k4
  
  )

dfmap <- listToMap(listmap) # make new sampleMap object

rownames(dfmap) <- dfmap$primary

# get object list
object.list <- list(
  
  snrnaseq.k2.all = sce1.all, 
  
  snrnaseq.k3.all = sce2.all, 
  
  snrnaseq.k4.all = sce3.all, 
  
  bulk.rnaseq = rse.filter,
  
  bulk.rpkm.rnaseq = rse.rpkm,
  
  cell.sizes = df.rnascope.kdata %>% as.matrix(),
  
  sce.img = sce.img,
  
  bulk.pb.k2 = list.pb.k234[["k2"]],
  
  bulk.pb.k3 = list.pb.k234[["k3"]],
  
  bulk.pb.k4 = list.pb.k234[["k4"]]
  
  )

experiment.list <- ExperimentList(object.list)

rownames(coldata) <- coldata$sample.id

# get mae datasets
mae <- prepMultiAssay(
  
  ExperimentList = experiment.list, 
                      
  
  sampleMap = dfmap, colData = coldata
  
  )

# get final mae type object
mae.final <- MultiAssayExperiment(
  
  mae$experiments, 
                                  
  
  mae$colData, 
                                  
  
  mae$sampleMap
  
  )



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

# inspect, with basic summaries
names(mae.final)
colData(mae.final)

#-----
# save
#-----

# set new mae filename
new.mae.filename <- "mae_allsamples.rda"

mae.final.filepath <- file.path("outputs", "01_mae", new.mae.filename)

save(mae.final, file = mae.final.filepath)
