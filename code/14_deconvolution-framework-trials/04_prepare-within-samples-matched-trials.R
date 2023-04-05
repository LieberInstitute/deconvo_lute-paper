#!/usr/bin/env R

# Author: Sean Maden
#
# Prepares experiment metadata for deconvolution within matched samples data.

source("deconvo_method-paper/code/14_deconvolution-framework-trials/00_parameters_script-set-14.R")
sapply(libv, library, character.only = T)
rse <- get(load(rse.k2markers.filepath))
image.table <- get(load(halo.output.path)) %>% as.data.frame()
image.cells <- get(load(image.cells.path))
sce <- get(load(sce.markers.list.path))[["k2"]]
# get logcounts expression
rse <- logNormCounts(rse, assay.type = "counts")
sce <- logNormCounts(sce, assay.type = "counts")

# get image-based cell size references
list.image.cell.sizes <- get(load(image.cell.sizes.save.path))
# get k2 cell size estimates
area.table <- list.image.cell.sizes$Nucleus_Area
counts.table <- list.image.cell.sizes$AKT3_Copies
colnames(area.table)[1] <- colnames(counts.table)[1] <- "cell.type"
colnames(area.table)[2] <- "nucleus.area"
colnames(counts.table)[2] <- "akt3.counts"
area.k2 <- get_k2_area(area.table, "nucleus.area")
area.k2 <- c("glial" = area.k2[1,2], "neuron" = area.k2[2,2])
counts.k2 <- get_k2_area(counts.table, "akt3.counts")
counts.k2 <- c("glial" = counts.k2[1,2], "neuron" = counts.k2[2,2])

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
# save
save(complete.sample.id.vector, file = complete.sample.id.vector.path)

# get list of deconvolution experiment object sets
# cell.sizes.manual <- c("glial" = 3, "neuron" = 10)
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
  image.cells.id <- image.cells[image.cells[,"sample.id"]==sample.id,]
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
