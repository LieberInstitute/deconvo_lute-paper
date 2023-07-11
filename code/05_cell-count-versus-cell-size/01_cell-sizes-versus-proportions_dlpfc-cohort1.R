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

# cell type harmonization table
df.ct <- data.frame(cell_type = c("Inhib", "Other", "Astro", "Endo", "Excit", "Oligo", "OPC", "Micro"),
                    k2 = c("neuron", "NA", "glial", "glial", "neuron", "glial", "glial", "glial"),
                    k3 = c("Inhib", "NA", "glial", "glial", "Excit", "glial", "glial", "glial"),
                    k4 = c("Inhib", "NA", "non_oligo_glial", "non_oligo_glial", "Excit", "Oligo", 
                           "non_oligo_glial", "non_oligo_glial"))

append_k_columns <- function(df.input, df.ct, celltype.variable = "cell_type"){
  celltype.vector <- df.input[,celltype.variable]
  unique.cell.types <- unique(celltype.vector)
  df.ct.filter <- df.ct[df.ct$cell_type %in% unique.cell.types,]
  df.ct.new <- do.call(rbind, lapply(celltype.vector, function(cell.type.id){
    df.ct.filter[df.ct.filter[,"cell_type"]==cell.type.id,,drop = F]
  }))
  cbind(df.input, df.ct.new)
}

# get cell size table
df.cellsize <- do.call(rbind, lapply(sample.id.vector, function(sample.id){
  mae.iter <- mae[,mae.cd$sample.id==sample.id,]
  
  # gets cell size
  img.data <- mae.iter[[3]] %>% assays() %>% unlist() %>% t() %>% as.data.frame()
  img.data$cell_type <- mae.iter[[3]]$cell_type
  df.iter <- aggregate(img.data[,"Nucleus_Area"], list("cell_type" = img.data$cell_type), FUN = "median")
  colnames(df.iter)[2] <- "median_nucleus_area"
  #df.iter$sample.id <- sample.id
  
  # gets cell amount by type, from sce
  stat <- sumCountsAcrossCells(mae.iter[[2]], ids = colData(mae.iter[[2]])[,"cellType_broad_hc"])
  stat.med <- colMedians(assays(stat)[["sum"]])
  df.stat <- data.frame(median_counts_snrnaseq = stat.med, cell_type = colnames(stat))
  # gets total cell amount by type, from sce
  stat2 <- sumCountsAcrossCells(mae.iter[[2]], ids = rep("all", ncol(mae.iter[[2]])))
  stat.med2 <- colMedians(assays(stat)[["sum"]])
  df.stat.all <- data.frame(median_counts_snrnaseq = stat.med2, cell_type = colnames(stat2))
  df.stat.final <- rbind(df.stat, df.stat.all)
  #df.stat.final$sample.id <- sample.id
  
  df.return <- cbind(df.stat.final, df.iter)
  df.return$sample.id <- sample.id
  df.return
}))


