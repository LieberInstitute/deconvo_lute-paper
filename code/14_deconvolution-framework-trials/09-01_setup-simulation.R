
# setup simulation
# perform deconvolution
# uses matched pseudobulk and signature matrices within samples
# includes null cell sizes adjustment

source("deconvo_method-paper/code/14_deconvolution-framework-trials/00_parameters_script-set-14.R")
sapply(libv, library, character.only = T)



sce <- get(load(sce.path))

# get logcounts expression
rse <- logNormCounts(rse, assay.type = "counts")
sce <- logNormCounts(sce, assay.type = "counts")