sce.name <- "sce_DLPFC.Rdata"
sce.path <- file.path("DLPFC_snRNAseq/processed-data/sce", sce.name)
sce.prepared.path <- file.path("deconvo_method-paper", "outputs", 
                               "15_sample-wise-signature-matrix-simulations",
                               "sce-prepared_dlpfc-ro1-train.rda")
sce <- get(load(sce.path))

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

sce <- scuttle::logNormCounts(sce)
markers.by.group <- markers_by_group(sce, 
                                     group.variable = "Sample", 
                                     celltype.variable = "k2", 
                                     assay.name = "logcounts", 
                                     markers.per.type = 500, 
                                     typemarker.algorithm = "meanratios",
                                     return.type = "list",
                                     verbose = TRUE)
markers.filtered <- filter_group_markers(markers.by.group,
                                         minimum.group.overlap.rate = 0.5)
# save overlap markers info
markers.filtered.path <- file.path("deconvo_method-paper/outputs/overlap-markers-filtered-k2_train.rda")
save(markers.filtered, file = markers.filtered.path)

# rank markers overall
markers.overall <- meanratiosParam(sce, "logcounts", "k2", 
                                   500, return.info = TRUE) %>% 
  typemarkers()
markers.overall.table <- markers.overall$result.info
filter.overall <- markers.overall.table$gene %in% 
  markers.filtered$marker.list.final
markers.overall.filtered <- markers.overall.table[filter.overall,]
markers.per.type.final <- 20
unique.cell.types <- markers.overall.filtered$cellType.target %>% unique()
markers.vector.overall.top <- lapply(unique.cell.types, 
                                     function(unique.type.id){
                                       markers.table.final <- markers.overall.filtered %>% 
                                         dplyr::filter(cellType.target == unique.type.id) %>% 
                                         dplyr::arrange(rank_ratio) %>% 
                                         dplyr::top_n(n = markers.per.type.final)
                                       markers.table.final$gene
                                     }) %>% unlist()

# save final markers overall
markers.vector.overall.top.path <- file.path("deconvo_method-paper/outputs/overall-markers-filtered-k2_train.rda")
save(markers.vector.overall.top, file = markers.vector.overall.top.path)

# save final markers se
filter.sce <- rownames(sce) %in% markers.vector.overall.top
sce <- sce[filter.sce,]
se.new <- SummarizedExperiment(assays = list(counts = as.matrix(counts(sce))),
                               colData = DataFrame(colData(sce)),
                               metadata = metadata(sce))
# save
se.markers.name <- "se-k2_markers-overlap-concordance-filter-final_train.rda"
se.markers.path <- here("deconvo_method-paper", "outputs", se.markers.name)
save(se.new, file = se.markers.path)

