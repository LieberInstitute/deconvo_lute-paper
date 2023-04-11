
# do overlap analyses of marker sets

script.path <- file.path("deconvo_method-paper", "code", 
                         "18_updated-marker-sets-analyses",
                         "01_load-marker-sets.R")
source(script.path)

# k2 marker types from sce data
se.markers.name <- "se-markers-overlap-filters_dlpfc-ro1-train.rda"
se.markers.path <- here("deconvo_method-paper", "outputs", 
                        "15_k2-simulations_within-sample-matched", se.markers.name)
se <- get(load(se.markers.path))
# get marker types
se <- se[rownames(se) %in% marker.list2$k2,]
sce <- SingleCellExperiment(se)
counts(sce) <- assays(se)[["counts"]]
sce <- scuttle::logNormCounts(sce)
z <- signature_matrix_from_sce(sce, "k2", "mean", "logcounts")


reference.list <- get(load(
  "~/GitHub/deconvo_method-paper/outputs/15_k2-simulations_within-sample-matched/list-references-by-overlap-rate_dlpfc-ro1-train.rda"))

ref.k2 <- reference.list[[1]][[1]]
length(intersect(rownames(ref.k2), marker.list2$k2))