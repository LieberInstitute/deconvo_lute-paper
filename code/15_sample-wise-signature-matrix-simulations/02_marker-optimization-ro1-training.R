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
min.overlap.rate <- NULL
markers.per.type.final <- 20
markers.per.type.iter.start <- 5000
# markers.by.type.iter <- seq(20, 1000, 20)

# cell.types.vector <- sce[["k2"]]
# unique.cell.types <- unique(cell.types.vector)

# get markers by overlaps
markers.by.batch <- markers_by_batch(sce = sce, "Sample", "k2", "logcounts", 
                                     markers.per.type.iter.start)
list.markers.final <- get.overlapping.markers(markers.by.batch,
                                              min.overlap.rate)
save(list.markers.final, file = list.markers.final.path)

list.markers.final <- get(load(here("deconvo_method-paper", "outputs", 
                                    "list-markers-by-slide-overlaps_dlpfc-ro1-train.rda")))

# filter duplicates
markers.vector <- lapply(list.markers.final, function(markers){markers[,1]}) %>% unlist()
markers.duplicated <- markers.vector[duplicated(markers.vector)]
list.markers.filtered <- lapply(list.markers.final, function(markers){
  markers[!markers[,1] %in% markers.duplicated,]
})

get_markers_overlap_filtered <- function(markers.list, min.overlap.rate = 80){
  # filter duplicates
  markers.vector <- lapply(list.markers.final, function(markers){markers[,1]}) %>% unlist()
  markers.duplicated <- markers.vector[duplicated(markers.vector)]
  list.markers.filtered <- lapply(list.markers.final, function(markers){
    markers[!markers[,1] %in% markers.duplicated,]
  })
  # get signature matrices for differing overlap frequencies
  markers.overlap.filtered <- lapply(names(list.markers.filtered), function(type){
    markers.table <- list.markers.filtered[[type]]
    overlap.filter <- markers.table$rate.overlap >= min.overlap.rate
    markers.table[overlap.filter,1]
  }) %>% unlist()
  # check duplicates
  message("found ", length(markers.overlap.filtered), " total markers.")
  message("found ", length(unique(markers.overlap.filtered)), " unique markers.")
  return(markers.overlap.filtered)
}

# get reference matrices
markers.list <- list.markers.filtered
overlap.rate.vector <- c(50, 60, 80, 90, 95)
list.reference.sets <- lapply(overlap.rate.vector, function(rate){
  get_markers_overlap_filtered(markers.list, rate)
})
reference.set.list <- lapply(list.reference.sets, function(set){
  sce.filter <- rownames(sce) %in% set
  sce.filtered <- sce[sce.filter,]
  message("getting the signature matrix for ", length(set), " markers.")
  signature_matrix_from_sce(sce.filtered,
                            cell.type.variable = "k2", 
                            summary.method = "mean", 
                            assay.name = "counts")
})

# get kappa results
# values by reference
kappa.list <- lapply(reference.set.list, function(reference){kappa(reference)})
# do resampling
num.rep <- 10
num.sample <- 2
kappa.list.sample2 <- lapply(reference.set.list, function(reference){
    sapply(seq(num.rep), 
           function(index){kappa(reference[sample(nrow(reference), 2),])})
})
median.resample <- sapply(kappa.list.sample2, median)
# values by assay matrix
kappa.list.sce <- lapply(reference.set.list, function(reference){
  sce[rownames(reference),] %>% counts() %>% kappa()
})

# append results
reference.list <- lapply(seq(length(reference.set.list)), function(index){
  list(reference = reference.set.list[[index]], 
       kappa.reference = kappa.list[[index]],
       kappa.sce = kappa.list.sce[[index]],
       kappa.resample2 = kappa.list.sample2[[index]],
       kappa.median.resample = median.resample[[index]])
})
# save results
names(reference.set.list) <- paste0("overlap_rate_", overlap.rate.vector)
reference.list.name <- "list-references-by-overlap-rate_dlpfc-ro1-train.rda"
reference.list.path <- here(save.path, reference.list.name)
save(reference.list, file = reference.list.path)
