#!/usr/bin/env R

# Author: Sean Maden
#
# Prepare cell quantities from halo image outputs.

source("deconvo_method-paper/code/14_deconvolution-framework-trials/00_parameters.R")
sapply(libv, library, character.only = T)
image.table <- get(load(halo.output.path))

# get cell amounts on brnum
image.brnum.vector <- unique(image.table[,"BrNum"])
image.cells <- do.call(rbind, lapply(image.brnum.vector, function(brnum){
  filter.brnum <- image.table[,"BrNum"]==brnum
  counts <- image.table[filter.brnum,"cell_type"] %>% table() 
  cell.types <- names(counts)
  counts <- counts %>% as.numeric()
  proportions <- counts %>% prop.table() %>% as.numeric()
  names(counts) <- names(proportions) <- cell.types
  names(counts) <- paste0(names(counts), ".count")
  names(proportions) <- paste0(names(proportions), ".proportion")
  c(counts, proportions, brnum)
})) %>% as.data.frame()
colnames(image.cells)[15] <- "BrNum"
for(index in seq(14)){
  image.cells[,index] <- image.cells[,index] %>% as.numeric()}
# get k2 values
column.names <- colnames(image.cells)
glial.id.vector <- c("Astro", "Oligo", "OPC", "Micro")
neuron.id.vector <- c("Inhib", "Excit")
which.glial.counts <- grepl(paste0(glial.id.vector, ".count", collapse = "|"), column.names)
which.neuron.counts <- grepl(paste0(neuron.id.vector, ".count", collapse = "|"), column.names)
image.cells$glial.count <- apply(image.cells[,which.glial.counts], 1, sum)
image.cells$neuron.count <- apply(image.cells[,which.neuron.counts], 1, sum)
image.cells$total.k2 <- apply(image.cells[,c("glial.count", "neuron.count")], 1, sum)
image.cells$glial.proportion <- image.cells$glial.count/image.cells$total.k2
image.cells$neuron.proportion <- image.cells$neuron.count/image.cells$total.k2
# save
save(image.cells, file = image.cells.path)
