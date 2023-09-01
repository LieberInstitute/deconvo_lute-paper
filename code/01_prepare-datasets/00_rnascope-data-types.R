#!/usr/bin/env R

#
# Gets RNAscope data formatted as SCE, and additional RNAscope data summaries.
#

libv <- c("here", "nlme", "ggplot2", "gridExtra", "dplyr", "ggforce", "MultiAssayExperiment", "SingleCellExperiment")
sapply(libv, library, character.only = TRUE)

#--------------------
# load prepped assays
#--------------------
# load rnascope image data
load("~/GitHub/deconvo_method-paper/outputs/01_prepare-datasets/halo-outputs_updated.Rdata")



#-----------------------------------
# rnascope: make sce with image data
#-----------------------------------
sample.id.halo <- "Sample"

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

#-------------------------------------
# absolute cell size filter on rowdata
#-------------------------------------
max.nucleus.size <- 75
filter.sce.size <- assays(sce.img)[["Nucleus_Area"]] <= max.nucleus.size
sce.img <- sce.img[,filter.sce.size]

#-----------------------------
# cell types filter on coldata
#-----------------------------
#cell.types.keep <- c("Excit", "Inhib", "Oligo")
#filter.sce.cd <- colData(sce.img)$cell_type %in% cell.types.keep
#sce.img <- sce.img[,filter.sce.cd]

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

#---------------
# save new files
#---------------
# sn.img
sce.img.path <- file.path("outputs", "01_prepare-datasets", "sce-rnascope_cohort1.rda")
save(sce.img, file = sce.img.path)
# df.rnascope.kdata
df.rn.path <- file.path("outputs", "01_prepare-datasets", "df-rnascope-info_cohort1.rda")
save(df.rnascope.kdata, file = df.rn.path)
