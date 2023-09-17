#!/usr/bin/env R
# Author: Sean Maden

source("./deconvo_method-paper/cohort1/scripts/04_experiment/00_parameters-script-set-15.R")
sapply(libv, library, character.only = T)
sce <- get(load(sce.path))

# functions
get_markers_overlap_filtered <- function(markers.list, 
                                         min.group.overlap = NULL, 
                                         verbose = TRUE){
  
  markers.list = markers.overlaps
  
  if(verbose){message("found ", length(unlist(markers.list)), " total markers.")}
  # filter duplicates
  markers.no.duplicates <- lapply(markers.list, function(type.markers){
    duplicated.type.markers <- type.markers[duplicated(type.markers)]
    type.markers[!type.markers[,1] %in% duplicated.type.markers,]
  })
  num.markers.unique <- length(unlist(markers.no.duplicates))
  if(verbose){message("found ", num.markers.unique, " non-duplicated markers.")}
  
  # filter concordant
  markers.unique.bytype <- lapply(markers.no.duplicates, 
                                  function(type.markers){type.markers[,1]}) %>% 
    unlist()
  markers.non.concordant <- markers.unique.bytype[
    duplicated(markers.unique.bytype)]
  markers.concordant <- lapply(markers.no.duplicates, function(type.markers){
    type.markers[!type.markers[,1] %in% markers.non.concordant,]
  })
  num.markers.non.concordant <- length(unlist(markers.non.concordant))
  if(verbose){message("found ", num.markers.non.concordant, " non-concordant markers.")}
  
  # get overlap rate filtered markers
  markers.overlap.filtered <- lapply(markers.concordant, function(type.markers){
    overlap.filter <- seq(nrow(type.markers))
    if(!is(min.overlap.rate, "NULL")){
      overlap.filter <- type.markers$rate.overlap >= min.overlap.rate
    }
    type.markers[overlap.filter, 1]
  }) %>% unlist() %>% unique()
  num.markers.overlap.filtered <- length(unlist(markers.overlap.filtered))
  if(verbose){
    message("found ", num.markers.overlap.filtered, 
            " group-overlapping markers with filter: ", min.group.overlap)}
  
  # return.list
  metadata.list <- list(num.markers.input.total = length(unlist(markers.list)),
                        num.markers.non.concordant = length(markers.non.concordant),
                        num.markers.overlap.filtered = length(markers.overlap.filtered),
                        min.group.overlap = min.group.overlap)
  return.list <- list(marker.list.final = markers.overlap.filtered,
                      filter.metadata = metadata.list)
  return(return.list)
}

# subset training donors
donor.variable <- "BrNum"
donor.id.train <- c("Br2720", "Br6471", "Br8492", "Br2743", "Br3942", "Br6423", "Br8325")
filter.train <- sce[[donor.variable]] %in% donor.id.train
sce <- sce[,filter.train]
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
markers.per.type.iter.start <- 10000
# markers.by.type.iter <- seq(20, 1000, 20)

# cell.types.vector <- sce[["k2"]]
# unique.cell.types <- unique(cell.types.vector)

# get markers by overlaps
markers.by.batch <- markers_by_batch(sce = sce, "Sample", 
                                     "k2", "logcounts", 
                                     markers.per.type.iter.start)

# save
markers.by.batch.path <- here("deconvo_method-paper", "cohort1", "data", 
                              "outputs", "04_experiment", 
                              "list-markers-by-batch_dlpfc-ro1-train.rda")
save(markers.by.batch, file = markers.by.batch.path)

markers.overlaps <- get.overlapping.markers(markers.by.batch, NULL)
markers.overlaps.path <- here("deconvo_method-paper", "cohort1", "data", 
                              "outputs", "04_experiment",
                              "list-markers-overlaps_dlpfc-ro1-train.rda")
save(markers.overlaps, file = markers.overlaps.path)

# apply additional filters
markers.filtered <- get_markers_overlap_filtered(markers.overlaps)
# save
markers.filtered.name <- "list-markers-overlap-filters_dlpfc-ro1-train.rda"
markers.filtered.path <- here("deconvo_method-paper", "cohort1", "data", 
                              "outputs", "04_experiment", markers.filtered.name)
save(markers.filtered, file = markers.filtered.name)

# get markers sce
filter.sce <- rownames(sce) %in% markers.filtered[[1]]
sce <- sce[filter.sce,]
dim(sce)
# append marker annotations
metadata(sce)$marker.data <- list(markers.overlaps = markers.overlaps,
                                  markers.filtered = markers.filtered)
# instantiate new sce object
se.new <- SummarizedExperiment(assays = list(counts = as.matrix(counts(sce))),
                                colData = DataFrame(colData(sce)),
                                metadata = metadata(sce))
# save
se.markers.name <- "se-markers-overlap-filters_dlpfc-ro1-train.rda"
se.markers.path <- here("deconvo_method-paper", "cohort1", "data", "outputs", 
                        "04_experiment", se.markers.name)
save(se.new, file = se.markers.path)

# new sce object
sce.new <- SingleCellExperiment(assays = list(counts = as.matrix(counts(sce))),
                                colData = DataFrame(colData(se.new)),
                                metadata = metadata(se.new))
# save
sce.markers.name <- "sce-markers-overlap-filters_dlpfc-ro1-train.rda"
sce.markers.path <- here("deconvo_method-paper", "cohort1", "data", "outputs", 
                         "04_experiment", sce.markers.name)
save(sce.new, file = sce.markers.path)

# get markers references
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

# save filtered sce object(s)
list.sce.data <- list(sce.filtered = sce.filtered, markers.data = markers.data)
list.sce.name <- ""
list.sce.path <- ""
save(list.sce.data, file = list.sce.path)

