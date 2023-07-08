#!/usr/bin/env R

# Author: Sean Maden
#
# Parameters for HALO image analyses.
#

# source("deconvo_method-paper/code/_halo-outputs-analysis/00_parameters.R")
# sapply(libv, library, character.only = T)


libv <- c("here", "nlme", "ggplot2", "gridExtra", "dplyr", "ggforce")
sapply(libv, library, character.only = TRUE)

# set variables
sample.id.label <- levels.variable <- "Sample"
cell.area.variable <- "Nucleus_Area"
gene.marker.label <- "AKT3_Copies"
normalized.area.variable <- "log10_nucleus-area"
normalized.marker.variable <- "log10_akt3-copies"

# set the halo data path
halo.output.file.name <- "halo_all.Rdata"
halo.output.path <- here("Human_DLPFC_Deconvolution", "processed-data", "03_HALO", halo.output.file.name)
halo.output.table <- get(load(halo.output.path))
halo.output.table <- halo.output.table %>% as.data.frame()

# normalize marker counts
# helper functions
normalization1 <- function(variable){log10(variable)}
# variable vectors
marker.vector <- halo.output.table[,gene.marker.label]
area.vector <- halo.output.table[,cell.area.variable]
# get log10 normalizations
halo.output.table[,normalized.area.variable] <- area.vector %>% normalization1() %>% as.numeric()
halo.output.table[,normalized.marker.variable] <- marker.vector %>% normalization1() %>% as.numeric()

# resave
# get save dpath
proj.dname <- "deconvo_method-paper"
code.dname <- "01_prepare-datasets"
save.path <- file.path(proj.dname, "outputs", code.dname)
output.updated.filename <- "halo-outputs_updated.Rdata"
output.updated.path <- here(save.path, output.updated.filename)
save(halo.output.table, file = output.updated.path)
