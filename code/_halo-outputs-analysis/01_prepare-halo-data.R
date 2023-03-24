#!/usr/bin/env R

# Author: Sean Maden
#
# Parameters for HALO image analyses.
#

source("deconvo_method-paper/code/_halo-outputs-analysis/00_parameters.R")
sapply(libv, library, character.only = T)
halo.output.table <- get(load(halo.output.path))
halo.output.table <- halo.output.table %>% as.data.frame()

# normalize marker counts
marker.vector <- halo.output.table[,gene.marker.label]
area.vector <- halo.output.table[,cell.area.variable]
# get log10 normalizations
halo.output.table[,normalized.area.variable] <- area.vector %>% normalization1() %>% as.numeric()
halo.output.table[,normalized.marker.variable] <- marker.vector %>% normalization1() %>% as.numeric()
# resave
save(halo.output.table, file = output.updated.path)
