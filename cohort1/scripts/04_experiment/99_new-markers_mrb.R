sce <- get(load(sce.path))
markers.train <- get(load(markers.train.path))

# key params
min.group.overlap.rate <- 0.5
markers.per.type.filter <- 500
markers.per.type.final <- 20
# save paths
markers.bygroup.info.path <- ""
markers.final.path <- ""
se.final.path <- ""

# prep sce
# get logcounts
sce <- scuttle::logNormCounts(sce)
# get celltype variable
sce[["k2"]] <- ifelse(grepl("Excit|Inhib", sce[["cellType"]]),"neuron",
                      ifelse(sce[["cellType"]] %in% c("Oligo", "OPC", "Micro", "Astro"), "glial", "other"))
# remove "other"
filter.sce <- !sce[["k2"]] == "other"
sce <- sce[,filter.sce]

# get markers by group
markers.by.group <- markers_by_group(sce, 
                 group.variable = "donor", 
                 celltype.variable = "k2", 
                 assay.name = "logcounts", 
                 markers.per.type = markers.per.type.filter, 
                 typemarker.algorithm = "meanratios",
                 return.type = "list",
                 verbose = FALSE)
markers.filtered <- filter_group_markers(markers.by.group,
                                         minimum.group.overlap.rate = 
                                           min.group.overlap.rate)

# rank markers overall
markers.overall <- meanratiosParam(sce, "logcounts", "k2", 
                                   markers.per.type.filter, 
                                   return.info = TRUE) %>% 
  typemarkers()
markers.overall.table <- markers.overall$result.info
filter.overall <- markers.overall.table$gene %in% 
  markers.filtered$marker.list.final
markers.overall.filtered <- markers.overall.table[filter.overall,]
unique.cell.types <- markers.overall.filtered$cellType.target %>% unique()
markers.vector.overall.top <- lapply(unique.cell.types, 
                                     function(unique.type.id){
  markers.table.final <- markers.overall.filtered %>% 
    dplyr::filter(cellType.target == unique.type.id) %>% 
    dplyr::arrange(rank_ratio) %>% 
    dplyr::top_n(n = markers.per.type.final)
  markers.table.final$gene
}) %>% unlist()

marker.overall.filtered.path <- "k2-markers_overall-filtered_mrb"
save(markers.vector.overall.top, file = marker.overall.filtered.path)

# compare markers
length(intersect(rownames(se.new), markers.vector.overall.top)) # [1] 349