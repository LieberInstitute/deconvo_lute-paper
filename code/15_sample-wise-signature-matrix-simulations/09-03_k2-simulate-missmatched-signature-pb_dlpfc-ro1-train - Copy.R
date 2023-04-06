
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
sce <- sce[,filter.train]
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
markers.vector <- rownames(reference.list[[3]][[1]])
sce <- sce[markers.vector,]

# experiment parameters
sample.id.variable <- "Sample"
celltype.variable <- "k2"
assay.name <- "counts"
s <- c("glial" = 3, "neuron" = 10)
deconvolution.algorithm <- "nnls"
summary.method  <- "mean"

experiment <- deconvolution.experiment.permute.groups(
  sce = sce, s = s,
  sample.id.variable = sample.id.variable, 
  celltype.variable = celltype.variable, 
  assay.name = assay.name,
  deconvolution.algorithm = deconvolution.algorithm, 
  summary.method = summary.method, 
  verbose = TRUE)

