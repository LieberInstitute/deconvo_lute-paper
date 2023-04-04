# do marker gene optimizations

source("deconvo_method-paper/code/15_sample-wise-signature-matrix-simulations/00_parameters-script-set-15.R")
sapply(libv, library, character.only = T)
sce <- get(load(sce.path))

# get k2 labels and filter
cell.type.vector <- sce[["cellType_broad_hc"]]
k2.cell.type.vector <- ifelse(cell.type.vector %in% c("Excit", "Inhib"), "neuron",
                              ifelse(cell.type.vector %in% c("Oligo", "OPC", "Astro", "Micro"), 
                                     "glial", "other"))
sce[["k2"]] <- k2.cell.type.vector
filter.other <- !sce[["k2"]]=="other"
sce <- sce[,filter.other]

# set params for overlaps
min.overlap.rate <- 0.8
markers.per.type.final <- 20
markers.per.type.iter.start <- 1000
# markers.by.type.iter <- seq(20, 1000, 20)

# cell.types.vector <- sce[["k2"]]
# unique.cell.types <- unique(cell.types.vector)

# get markers by overlaps
markers.by.batch <- markers_by_batch(sce = sce, "Sample", "k2", "logcounts", 
                                     markers.per.type.iter.start)
list.markers.final <- get.overlapping.markers(markers.by.batch,
                                              min.overlap.rate)
list.markers.final






if(nrow(markers.final))




