
# setup simulation
# perform deconvolution
# uses matched pseudobulk and signature matrices within samples
# includes null cell sizes adjustment

source("deconvo_method-paper/code/15_sample-wise-signature-matrix-simulations/00_parameters-script-set-15.R")
sapply(libv, library, character.only = T)

# get save dpath
sce.path <- here("deconvo_method-paper", "outputs", "09_manuscript", "sce-mrb_dlpfc.rda")
sce <- get(load(sce.path))

# get k2 labels and filter
cell.type.vector <- sce[["cellType"]]
k2.cell.type.vector <- ifelse(grepl("Excit|Inhib", cell.type.vector), "neuron",
                              ifelse(cell.type.vector %in% c("Oligo", "OPC", "Astro", "Micro"), 
                                     "glial", "other"))
sce[["k2"]] <- k2.cell.type.vector
filter.other <- !sce[["k2"]]=="other"
sce <- sce[,filter.other]

# get concordant, overlapping markers
#save.path <- here("deconvo_method-paper", "outputs")
#reference.list.name <- "list-references-by-overlap-rate_dlpfc-ro1-train.rda"
reference.list.path <- here("deconvo_method-paper", "outputs", 
                            "15_k2-simulations_within-sample-matched", 
                            "list-references-by-overlap-rate_dlpfc-ro1-train.rda")
reference.list <- get(load(reference.list.path))
# get markers for this cell type
markers.vector <- rownames(reference.list[[3]][[1]])
length(markers.vector)

# experiment parameters
assay.name <- "counts"
celltype.variable <- "k2"
sample.id.variable <- "donor"
S.pb <- c("glial" = 3, "neuron" = 10)
deconvolution.algorithm <- "nnls"
experiment <- deconvolution.experiment(y = NULL, sce = sce, s = S.pb, 
                                       assay.name = assay.name,
                                       sample.id.variable = sample.id.variable, 
                                       celltype.variable = celltype.variable, 
                                       markers.vector = markers.vector)

# view plots
experiment$plots.list$neuron$proportions.scatterplot
experiment$plots.list$neuron$abs.error.barplot
experiment$plots.list$neuron$abs.error.jitterbox

# save
experiment.name <- "experiment_within-slide-matched_dlpfc-mrb.rda"
experiment.path <- here("deconvo_method-paper", "outputs", experiment.name)
save(experiment, file = experiment.path)
