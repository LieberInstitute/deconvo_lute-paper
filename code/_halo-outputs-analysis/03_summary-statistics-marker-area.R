#!/usr/bin/env R

# Author: Sean Maden
#

source("deconvo_method-paper/code/_halo-outputs-analysis/00_parameters.R")
sapply(libv, library, character.only = T)
halo.outputs.table <- get(load(halo.output.path))
halo.outputs.table <- halo.outputs.table %>% as.data.frame()
cell.area.vector <- halo.outputs.table[,cell.area.variable]
gene.marker.vector <- halo.outputs.table[,gene.marker.label]
sample.id.vector <- halo.outputs.table[,sample.id.label]

summary.matrix <- do.call(rbind, lapply(unique(sample.id.vector), function(sample.id){
  # get filtered data
  filter <- halo.outputs.table[,sample.id.label]==sample.id
  table.filtered <- halo.outputs.table[filter,]
  # get filtered variable data
  cell.type.list <- list(cell.type = table.filtered[, cell.type.label])
  cell.area.filtered <- table.filtered[, cell.area.variable]
  marker.counts.filtered <- table.filtered[, gene.marker.label]
  # get matrix summary statistics
  do.call(rbind, lapply(summary.terms, function(term){
    matrix.area <- aggregate(cell.area.filtered, by = cell.type.list, FUN = term)
    matrix.marker <- aggregate(marker.counts.filtered, by = cell.type.list, FUN = term)
    colnames(matrix.area)[2] <- colnames(matrix.marker)[2] <- "value"
    matrix.area$metric <- "nucleus_area"
    matrix.marker$metric <- "akt3_copies"
    matrix.area$summary.term <- matrix.marker$summary.term <- term
    matrix.marker$sample.id <- matrix.area$sample.id <- sample.id
    rbind(matrix.area, matrix.marker)
  }))
}))
