
# setup simulation
# perform deconvolution
# uses missmatched pseudobulk and signature matrices within samples
# includes null cell sizes adjustment

source("deconvo_method-paper/code/15_sample-wise-signature-matrix-simulations/00_parameters-script-set-15.R")
sapply(libv, library, character.only = T)

# load sce
sce.path <- here("DLPFC_snRNAseq/processed-data/sce", "sce_DLPFC.Rdata")
sce <- get(load(sce.path))
# subset training donors
donor.variable <- "BrNum"
donor.id.train <- c("Br2720", "Br6471", "Br8492", "Br2743", "Br3942", "Br6423", "Br8325")
filter.train <- sce[[donor.variable]] %in% donor.id.train
sce.train <- sce[,filter.train]
# get k2 labels and filter
cell.type.vector <- sce[["cellType_broad_hc"]]
k2.cell.type.vector <- ifelse(cell.type.vector %in% c("Excit", "Inhib"), "neuron",
                              ifelse(cell.type.vector %in% c("Oligo", "OPC", "Astro", "Micro"), 
                                     "glial", "other"))
sce[["k2"]] <- k2.cell.type.vector
filter.other <- !sce[["k2"]]=="other"
sce <- sce[,filter.other]

# filter markers
#reference.list.path <- here("deconvo_method-paper", "outputs", 
#                            "15_k2-simulations_within-sample-matched", 
#                            "list-references-by-overlap-rate_dlpfc-ro1-train.rda")
#reference.list <- get(load(reference.list.path))
#markers.vector <- rownames(reference.list[[3]][[1]])
#sce <- sce[markers.vector,]
save.path <- here("deconvo_method-paper", "outputs")
reference.list.name <- "list-references-by-overlap-rate_dlpfc-ro1-train.rda"
reference.list.path <- here(save.path, reference.list.name)
reference.list <- get(load(reference.list.path))

# experiment parameters
sample.id.variable <- "donor"
celltype.variable <- "k2"
assay.name <- "counts"
s <- c("glial" = 3, "neuron" = 10)
deconvolution.algorithm <- "nnls"
summary.method  <- "mean"

# get sample/group ids
group.id.vector <- sce[[sample.id.variable]]
unique.group.id.vector <- group.id.vector %>% unique() %>% as.character()
results.table.list <- lapply(unique.group.id.vector, function(group.id.z){
  
  message("working on group id: ", group.id.z)
  message("getting main signature matrix...")
  filter.sce <- sce[[sample.id.variable]]==group.id.z
  sce.group.z <- sce[,filter.sce]
  z <- signature_matrix_from_sce(sce.group.z,
                                 celltype.variable = celltype.variable,
                                 summary.method = summary.method,
                                 assay.name = assay.name)
  
  message("getting pseudobulks...")
  pb.filter <- !unique.group.id.vector==group.id.z
  unique.group.id.pb <- unique.group.id.vector[pb.filter] %>% unique()
  ypb.list <- lapply(unique.group.id.pb, function(group.id){
    message("working on pseudobulk for sample id ", group.id, "...")
    filter.group <- sce[[sample.id.variable]]==group.id
    ypb_from_sce(sce = sce[,filter.group], assay.name = assay.name, 
                 celltype.variable = celltype.variable, S = s)
  })
  ypb.table <- do.call(cbind, ypb.list) %>% as.data.frame()
  colnames(ypb.table) <- unique.group.id.pb
  
  message("getting experiment series...")
  experiment <- deconvolution.experiment(sce = sce, 
                                         y = ypb.table, s = s, z = z, 
                                         assay.name = assay.name, 
                                         sample.id.variable = sample.id.variable, 
                                         experiment.labels = group.id.z,
                                         deconvolution.algorithm = deconvolution.algorithm,
                                         celltype.variable = celltype.variable)
  
  message("finished with group id, returning results table...")
  # get results table
  results.table.iteration <- experiment$results.table
  results.table.iteration$group.id.signature <- group.id.z
  results.table.iteration
})
results.table <- do.call(rbind, results.table.list) %>% as.data.frame()

# get plots
# all results together
lgg <- deconvolution.results.plots(results.table)
# just miss-matched simulations
filter.results <- results.table$sample.id==results.table$group.id.signature
results.filtered <- results.table[!filter.results,]
results.filtered$sample.id <- results.filtered$group.id.signature
lgg.missmatch <- deconvolution.results.plots(results.filtered)
# just matched simulations
results.filtered <- results.table[filter.results,]
results.filtered$sample.id <- results.filtered$group.id.signature
lgg.match <- deconvolution.results.plots(results.filtered)

# print plots
# missmatched
title.string <- "Neuron, missmatched Y, Z"
lgg.missmatch$neuron$proportions.scatterplot + ggtitle(title.string)
lgg.missmatch$neuron$abs.error.jitterbox + facet_wrap(~sample.id) + 
  ggtitle(title.string)
# matched
title.string <- "Neuron, matched Y, Z"
lgg.match$neuron$proportions.scatterplot + ggtitle(title.string)
lgg.match$neuron$abs.error.barplot + ggtitle(title.string)
lgg.match$neuron$abs.error.jitterbox + facet_wrap(~sample.id) + 
  ggtitle(title.string)





