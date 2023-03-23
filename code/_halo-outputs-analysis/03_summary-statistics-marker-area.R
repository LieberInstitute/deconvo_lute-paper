#!/usr/bin/env R

# Author: Sean Maden
#

source("deconvo_method-paper/code/_halo-outputs-analysis/00_parameters.R")
sapply(libv, library, character.only = T)
halo.outputs.table <- get(load(halo.output.path))
halo.outputs.table <- halo.outputs.table %>% as.data.frame()

# check equality of variance between groups
cell.area.vector <- halo.outputs.table[,cell.area.variable]
gene.marker.vector <- halo.outputs.table[,gene.marker.label]
sample.id.vector <- halo.outputs.table[,sample.id.label]


summary_term_list(table = halo.outputs.table, 
                  sample.id.vector = sample.id.label, 
                  summary.variable.label = ,
                  summary.terms = c("variance", "mean", "max", "min"))

marker.variance.matrix <- do.call(cbind, lapply(sample.id.vector, function(sample.id){
  filter <- halo.outputs.table[,sample.id.label] == sample.id
  halo.outputs.table[filter, gene.marker.label] %>% as.numeric() %>% var()
}))
