#!/usr/bin/env R

# Author: Sean Maden
#
#

libv <- c("SingleCellExperiment", "SummarizedExperiment")
sapply(libv, library, character.only = T)

#----------
# load data
#----------
# load bulk data
rse.filename <- "rse_gene.Rdata"
load.path <- file.path("Human_DLPFC_Deconvolution", "processed-data",
                       "01_SPEAQeasy", "round2_v40_2022-07-06", "rse")
rse <- get(load(file.path(load.path, rse.filename)))

# get save directory path
save.path <- file.path("deconvo_method-paper", "outputs", 
                       "11_bulk-expression-analysis")

# load marker data
sce.markers.path <- file.path("deconvo_method-paper", "outputs", "09_manuscript", 
                              "list-scef_markers-k2-k3-k4_ro1-dlpfc.rda")
lsce <- get(load(sce.markers.path))
sce <- lsce[["k2"]]
rm(lsce)

#--------------------------------
# get marker gene bulk expression
#--------------------------------
# subset rse
bulk.gene.names <- rowData(rse)$Symbol
marker.genes.vector <- rownames(sce)
overlapping.markers <- intersect(bulk.gene.names, marker.genes.vector)
message("Found ", length(overlapping.markers), " overlapping markers.")
filter <- which(rowData(rse)$Symbol %in% overlapping.markers)
rsef <- rse[filter,]
dim(rsef)

# save bulk marker expr
rsef.filename <- "rsef_k2-marker-expr_ro1-dlpfc.rda"
save.path <- file.path(save.path, rsef.filename)
save(rsef, file = save.path)

#--------------------------------
# qc at markers versus background
#--------------------------------







