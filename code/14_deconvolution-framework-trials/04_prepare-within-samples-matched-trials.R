#!/usr/bin/env R

# Author: Sean Maden
#
# Prepares experiment metadata for deconvolution within matched samples data.

source("deconvo_method-paper/code/14_deconvolution-framework-trials/00_parameters.R")
sapply(libv, library, character.only = T)
rse <- get(load(rse.k2markers.filepath))
image.table <- get(load(halo.output.path))
image.cells <- get(load(image.cells.path))
sce <- get(load(sce.markers.list.path))[["k2"]]
sce <- logNormCounts(sce)

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
intersect(unique.sample.id.sce, unique.sample.id.image)
intersect(unique.sample.id.sce, unique.sample.id.rse)
intersect(unique.sample.id.rse, unique.sample.id.image)
complete.sample.id.vector <- intersect(unique.sample.id.sce, 
                                       intersect(unique.sample.id.image, 
                                                 unique.sample.id.rse))
# get list of deconvolution experiment object sets
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
  image.cells.id <- image.cells[image.cells[,"BrNum"]==brnum,]
  p.count.k2 <- c(image.cells.id["glial.count"], 
                  image.cells.id["neuron.count"]) %>% unlist()
  p.proportion.k2 <- c(image.cells.id["glial.proportion"], 
                       image.cells.id["neuron.proportion"]) %>% unlist()
  names(p.count.k2) <- names(p.proportion.k2) <- c("glial", "neuron")
  y <- assays(rse.sample)[[assay.name.rse]]
  return(list(sce = sce.sample, image = image.sample, rse = rse.sample,
              z = z.sample, p.count.k2 = p.count.k2, 
              p.proportion.k2 = p.proportion.k2, y = y))
})
names(lexperiment) <- complete.sample.id.vector
# save lexperiment
save(lexperiment, file = lexperiment.withinsample.path)
