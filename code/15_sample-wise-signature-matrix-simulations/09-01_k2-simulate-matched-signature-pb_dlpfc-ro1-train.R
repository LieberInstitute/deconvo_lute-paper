
# setup simulation
# perform deconvolution
# uses matched pseudobulk and signature matrices within samples
# includes null cell sizes adjustment

source("deconvo_method-paper/code/15_sample-wise-signature-matrix-simulations/00_parameters_script-set-15.R")
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

# get concordant, overlapping markers
#reference.list.path <- here("deconvo_method-paper", "outputs", 
#                            "15_k2-simulations_within-sample-matched", 
#                            "list-references-by-overlap-rate_dlpfc-ro1-train.rda")
save.path <- here("deconvo_method-paper", "outputs")
reference.list.name <- "list-references-by-overlap-rate_dlpfc-ro1-train.rda"
reference.list.path <- here(save.path, reference.list.name)
reference.list <- get(load(reference.list.path))
# get markers for this cell type
markers.vector <- rownames(reference.list[[3]][[1]])
length(markers.vector)

# experiment parameters
assay.name <- "counts"
celltype.variable <- "k2"
sample.id.variable <- "Sample"
S.pb <- c("glial" = 3, "neuron" = 10)
deconvolution.algorithm <- "nnls"
experiment <- deconvolution_experiment(y = NULL, sce = sce, s = S.pb, 
                                       assay.name = assay.name,
                                       sample.id.variable = sample.id.variable, 
                                       celltype.variable = celltype.variable, 
                                       markers.vector = markers.vector)

# view plots
experiment$plots.list$neuron$proportions.scatterplot
experiment$plots.list$neuron$abs.error.barplot
experiment$plots.list$neuron$abs.error.jitterbox

# save
experiment.name <- "experiment_within-slide-matched_dlpfc-ro1-train.rda"
experiment.path <- here("deconvo_method-paper", "outputs", experiment.name)
save(experiment, file = experiment.path)
