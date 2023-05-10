# dependencies
libv <- c("here", "lute", "dplyr", "ggplot2", "gridExtra", "SingleCellExperiment", 
          "SummarizedExperiment", "scran", "DeconvoBuddies", "UpSetR", "DelayedArray")
sapply(libv, library, character.only = TRUE)

# save path
save.path <- here("deconvo_method-paper", "outputs", "20_snrnaseq-bulk-matched_training")
save.path <- here("deconvo_method-paper", "outputs")

# get expression data
# snrnaseq
# bulk rnaseq

# get markers
# save markers

# get deconvolution results for matched samples
# without sfactor transformation
s <- c(1, 1)
# with some sfactor transformation
s <- c(3, 10)
