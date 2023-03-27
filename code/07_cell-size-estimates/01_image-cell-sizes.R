#!/usr/bin/env R

# Author: Sean Maden
#
# Get cell sizes from formatted image analysis outputs.

source("deconvo_method-paper/code/07_cell-size-estimates/00_parameters.R")
sapply(libv, library, character.only = T)
halo.output.table <- get(load(halo.output.path))
halo.output.table <- halo.output.table %>% as.data.frame()
# get cell size estimates
list.image.cell.sizes <- lapply(area.variables, function(area.variable){
  sizes <- cell_sizes(halo.output.table, area.variable = area.variable)
  colnames(sizes) <- c("cell_type", "size")
  sizes$area.variable <- area.variable
  sizes
})
names(list.image.cell.sizes) <- area.variables
# save
save(list.image.cell.sizes, file = image.cell.sizes.save.path)