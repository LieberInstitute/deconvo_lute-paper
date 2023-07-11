# dependencies
libv <- c("here", "lute", "dplyr", "ggplot2", "gridExtra", "SingleCellExperiment", 
          "SummarizedExperiment", "scran")
sapply(libv, library, character.only = TRUE)

# save path
prepped.data.path <- here("deconvo_method-paper", "outputs", "01_prepare-datasets")
save.path <- here("deconvo_method-paper", "outputs", "05_cell-count-versus-cell-size")

# multi assay experiment path
mae.filename <- ""
mae.path <- here(save.path, )