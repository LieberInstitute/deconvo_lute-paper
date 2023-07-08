
# setup simulation
# perform deconvolution
# uses missmatched pseudobulk and signature matrices within samples
# includes null cell sizes adjustment

source("deconvo_method-paper/code/15_sample-wise-signature-matrix-simulations/00_parameters-script-set-15.R")
sapply(libv, library, character.only = T)
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

# filter markers
reference.list.path <- here("deconvo_method-paper", "outputs", 
                            "15_k2-simulations_within-sample-matched", 
                            "list-references-by-overlap-rate_dlpfc-ro1-train.rda")
reference.list <- get(load(reference.list.path))
markers.vector <- rownames(reference.list[[3]][[1]])
sce <- sce[markers.vector,]

# experiment parameters
sample.id.variable <- "donor"
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
