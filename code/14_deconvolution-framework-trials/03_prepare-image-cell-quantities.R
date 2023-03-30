#!/usr/bin/env R

# Author: Sean Maden
#
# Prepare cell quantities from halo image outputs.

source("deconvo_method-paper/code/14_deconvolution-framework-trials/00_parameters.R")
sapply(libv, library, character.only = T)
image.table <- get(load(halo.output.path))

# get k2 labels
cell.type.vector <- image.table$cell_type
image.table$k2.type <- ifelse(cell.type.vector %in% c("Inhib", "Excit"), "neuron",
                              ifelse(cell.type.vector %in% 
                                       c("Oligo", "OPC", "Astro", "Micro"), 
                                     "glial", "other"))
# filter other
image.table <- image.table[!image.table$k2.type=="other",]

# get proportions on slide (combine circle/star)
image.table$sample.id <- paste0(image.table$BrNum, "_", toupper(image.table$Position))
unique.sample.id <- unique(image.table$sample.id)
image.prop <- do.call(rbind, lapply(unique.sample.id, function(sample.id){
  filter <- image.table$sample.id == sample.id
  table(image.table[filter,]$k2.type) %>% prop.table() %>% t()
})) %>% as.data.frame()
image.prop$sample.id <- unique.sample.id
# save
save(image.prop, file = image.cells.path)
