#!/usr/bin/env R

#
# Testing relationship between cell sizes and total cells.
#

# aggregation basis code:
# do.call(rbind, lapply(sample.id.vector, function(sample.id){}))

source("deconvo_method-paper/code/05_cell-count-versus-cell-size/00_parameters.R")
sapply(libv, library, character.only = T)

# load mae 
mae.filepath <- here("deconvo_method-paper", "outputs", "01_prepare-datasets", "mae_final.rda")
mae <- get(load(mae.filepath))
# get mae metadata for subsetting
mae.cd <- colData(mae.final)
sample.id.vector <- mae.cd$sample.id[complete.cases(mae.final)]

# get cell size table
df.cellsize <- do.call(rbind, lapply(sample.id.vector, function(sample.id){
  mae.iter <- mae[,mae.cd$sample.id==sample.id,]
  img.data <- mae.iter[[3]] %>% assays() %>% unlist() %>% t() %>% as.data.frame()
  img.data$cell_type <- mae.iter[[3]]$cell_type
  df.iter <- aggregate(img.data[,"Nucleus_Area"], list("cell_type" = img.data$cell_type), FUN = "median")
  colnames(df.iter)[2] <- "median_nucleus_area"
  df.iter$sample.id <- sample.id
  df.iter
}))


