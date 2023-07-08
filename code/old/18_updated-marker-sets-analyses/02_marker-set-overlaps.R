
# do overlap analyses of marker sets
# compare marker sets, with and without concordance/overlap filters
# k2 markers only

script.path <- file.path("deconvo_method-paper", "code", "18_updated-marker-sets-analyses",
                         "01_load-marker-sets.R")
source(script.path)
rm(se.filter)
rm(sutton.format)

# get top 40 concordant/overlapping
# markers overlapping
index.vector <- seq(length(marker.train.group))
marker.group.input <- lapply(index.vector, function(index){
  data.table <- marker.train.group[[index]]$result.info
  data.table$group.id <- names(marker.train.group)[index]
  return(data.table)
})
names(marker.group.input) <- names(marker.train.group)
markers.filter <- filter_group_markers(marker.group.input, 0.5)

filter.se <- rownames(se.train) %in% markers.filter$marker.list.final
se.train.filter <- se.train[filter.se,]
dim(se.train.filter)
sce <- SingleCellExperiment(se.train.filter)
assays(sce) <- list(counts = assays(se.train.filter)[["counts_adjusted"]])
sce <- scuttle::logNormCounts(sce)
k2.mrb.top40.filter <- meanratiosParam(sce = sce, assay.name = "logcounts", 
                                       celltype.variable = "k2", return.info = T,
                                       markers.per.type = 20) %>% typemarkers()
  

# marker overlap
k2.train.top40.nofilter <- metadata(se.train)$top.40.markers.results[[1]]
k2.mrb.top40.nofilter <- metadata(se.mrb)$top.40.markers.results[[1]]



k2.train.top40.filter <- intersect(k2.train.top40.nofilter, 
                                   marker.train.overlaps$marker.list.final)
k2.mrb.top40.filter <- intersect(k2.mrb.top40.nofilter, 
                                   marker.mrb.overlaps$marker.list.final)


# mrb k2 markers expression
se.markers.name <- "se-markers-overlap-filters_dlpfc-ro1-train.rda"
se.markers.path <- here("deconvo_method-paper", "outputs", 
                        "15_k2-simulations_within-sample-matched", se.markers.name)
se <- get(load(se.markers.path))

# ro1 k2 markers expression
# get marker types
se <- se[rownames(se) %in% marker.list2$k2,]
sce <- SingleCellExperiment(se)
counts(sce) <- assays(se)[["counts"]]
sce <- scuttle::logNormCounts(sce)
z <- signature_matrix_from_sce(sce, "k2", "mean", "logcounts")
