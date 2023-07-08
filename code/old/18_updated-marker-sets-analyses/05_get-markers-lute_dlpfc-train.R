
libv <- c("here", "dplyr", "lute", "biomaRt", "SummarizedExperiment", 
          "SingleCellExperiment", "scuttle")
sapply(libv, library, character.only = TRUE)

{
  # main run parameters
  markers.per.group.discover <- 500
  markers.per.group.final <- 40
  celltype.variable.original <- "cellType_broad_hc"
  celltype.variable.new <- "k2"
  group.id.variable <- "Sample"
  assay.name.markers <- "counts"
  typemarker.algorithm <- "meanratios"
  min.group.overlap.rate <- 0.5
  markers.bygroup.name <- "group-markers-k2_train.rda"
  markers.filtered.name <- "overlap-markers-filtered-k2_train.rda"
  
  # get sce data
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
  cell.type.vector <- sce[[celltype.variable.original]]
  k2.cell.type.vector <- ifelse(cell.type.vector %in% c("Excit", "Inhib"), "neuron",
                                ifelse(cell.type.vector %in% c("Oligo", "OPC", "Astro", "Micro"), 
                                       "glial", "other"))
  sce[[celltype.variable.new]] <- k2.cell.type.vector
  filter.other <- !sce[[celltype.variable.new]]=="other"
  sce <- sce[,filter.other]
  sce <- scuttle::logNormCounts(sce)
}


# get all marker candidates
{
  markers.by.group <- markers_by_group(sce, 
                                       group.variable = group.id.variable, 
                                       celltype.variable = celltype.variable.new, 
                                       assay.name = "logcounts", 
                                       markers.per.type = markers.per.group.discover, 
                                       typemarker.algorithm = typemarker.algorithm,
                                       return.type = "list",
                                       verbose = TRUE)
}
