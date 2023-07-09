#!/usr/bin/env R

# Author: Sean Maden
#
# Prepares experiment metadata for deconvolution within matched samples data.

source("deconvo_method-paper/code/03_matched-bulk-predictions/00_parameters.R")
sapply(libv, library, character.only = T)

# load
rse <- get(load(rse.k2markers.filepath))
image.table <- get(load(halo.output.path)) %>% as.data.frame()
sce <- get(load(sce.markers.list.path))[["k2"]]
halo.cellsize <- get(load(halo.cellsize.filepath))
halo.count <- get(load(halo.cellcount.filepath))
halo.prop <- get(load(halo.cellprop.filepath))

# get logcounts expression
rse <- logNormCounts(rse, assay.type = "counts")
sce <- logNormCounts(sce, assay.type = "counts")

# get matched datasets
# image
image.sample.id.vector <- standard_sample_id(image.table, "Position")
sce.sample.id.vector <- standard_sample_id(colData(sce), "Position")
rse.sample.id.vector <- standard_sample_id(colData(rse), "location")
# unique sample ids
unique.sample.id.sce <- unique(sce.sample.id.vector)
unique.sample.id.rse <- unique(rse.sample.id.vector)
unique.sample.id.image <- unique(image.sample.id.vector)
# overlapping sample ids
length(intersect(unique.sample.id.sce, unique.sample.id.image)) # 11
length(intersect(unique.sample.id.sce, unique.sample.id.rse)) # 11
length(intersect(unique.sample.id.rse, unique.sample.id.image)) # 19
complete.sample.id.vector <- intersect(unique.sample.id.sce, 
                                       intersect(unique.sample.id.image, 
                                                 unique.sample.id.rse))

#-------------------------------------------------
# get list of deconvolution experiment object sets
#-------------------------------------------------
lexperiment <- lapply(complete.sample.id.vector, function(sample.id){
  message(sample.id)
  filter.sce <- sce.sample.id.vector == sample.id
  filter.image <- image.sample.id.vector == sample.id
  filter.rse <- rse.sample.id.vector == sample.id
  # get sample data
  image.sample <- image.table[filter.image,]
  sce.sample <- sce[,filter.sce]
  rse.sample <- rse[,filter.rse]
  
  
  # get deconvolution objects
  z.sample <- signature_matrix_from_sce(sce.sample)
  
  # parse cell counts/proportions
  brnum <- gsub("_.*", "", sample.id)
  
  # filtered cellsizes
  halo.cellsize.sample <- lapply(halo.cellsize, function(sample.label.data){
    data.sample <- sample.label.data[["k2"]]
    data.sample$sample.id <- paste0(gsub("_.*", "", data.sample$sample.id), "_", 
                               toupper(gsub(".*_", "", data.sample$sample.id)))
    sample.id.fetch.format <- paste0(gsub("_.*", "", sample.id), "_",
                                     paste0(unlist(strsplit(gsub(".*_", "", sample.id),""))[1:3],collapse =""))
    data.sample[data.sample$sample.id==sample.id.fetch.format,]
  })
  
  halo.count.sample <- lapply(halo.count, function(sample.label.data){
    data.sample <- sample.label.data[["k2"]]
    data.sample$sample.id <- paste0(gsub("_.*", "", data.sample$sample.id), "_", 
                                    toupper(gsub(".*_", "", data.sample$sample.id)))
    sample.id.fetch.format <- paste0(gsub("_.*", "", sample.id), "_",
                                     paste0(unlist(strsplit(gsub(".*_", "", sample.id),""))[1:3],collapse =""))
    data.sample[data.sample$sample.id==sample.id.fetch.format,]
  })
  
  p.count.k2 <- c(image.cells.id["glial.count"], image.cells.id["neuron.count"]) %>% unlist()
  p.proportion.k2 <- c(image.cells.id["glial.proportion"], 
                       image.cells.id["neuron.proportion"]) %>% unlist()
  
  names(p.count.k2) <- names(p.proportion.k2) <- c("glial", "neuron")
  
  y <- assays(rse.sample)[[assay.name.rse]]
  colnames(y) <- rse.sample[["expt_condition"]]
  
  # set cell sizes
  # sizes.sample <- cell_size_sample(image.table, sample.id, image.sample.id.vector)
  list.sizes <- list(reference.area = area.k2, 
                     reference.counts = counts.k2,
                     null = c("glial" = 1, "neuron" = 1))
  # get total expression y
  y.total.expression <- colSums(y)
  # return
  return(list(sce = sce.sample, image = image.sample, rse = rse.sample,
              z = z.sample, p.count.k2 = p.count.k2, 
              p.proportion.k2 = p.proportion.k2, y = y,
              y.total.expression = y.total.expression,
              list.sizes = list.sizes))
})
names(lexperiment) <- complete.sample.id.vector
# save lexperiment
save(lexperiment, file = lexperiment.withinsample.path)
