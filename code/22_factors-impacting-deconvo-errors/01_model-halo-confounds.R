

#
# revisits ancova for halo confounds, with the following revisions:
# * subsample cells such that residuals are minimal
# * take subsets of cell sizes by size ranges such that residuals are minimal
# * do pca to compare cells, samples
#
#

# dependencies
libv <- c("here", "lute", "dplyr", "ggplot2", "gridExtra", "SingleCellExperiment", "scuttle",
          "SummarizedExperiment", "scran", "DeconvoBuddies", "UpSetR", "DelayedArray")
sapply(libv, library, character.only = TRUE)

# save path
save.path <- here("deconvo_method-paper", "outputs", "20_snrnaseq-bulk-matched_training")

#--------------
# get halo data
#--------------
# image reference
halo.path <- "Human_DLPFC_Deconvolution/processed-data/03_HALO/halo_all.Rdata"
halo.all <- get(load(halo.path))
