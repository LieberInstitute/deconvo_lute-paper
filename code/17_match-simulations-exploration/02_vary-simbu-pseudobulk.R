# simulations using slide-matched y,z, but varying the simbu-introduced 
# pseudobulk biases in reflection of true bulk datasets.

libv <- c("here", "dplyr", "lute", "SummarizedExperiment", 
          "SingleCellExperiment", "ggplot2", "ggcorrplot",
          "reshape2", "ggrepel", "gridExtra")
sapply(libv, library, character.only = TRUE)

save.path <- here("deconvo_method-paper", "outputs", 
                  "17_match-simulations-exploration")