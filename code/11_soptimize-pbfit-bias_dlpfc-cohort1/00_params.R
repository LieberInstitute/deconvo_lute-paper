#!/usr/bin/env R

# Author: Sean Maden

# dependencies
libv <- c("lute", "dplyr", "ggplot2", "gridExtra", "SingleCellExperiment", "SummarizedExperiment", "scran")
sapply(libv, library, character.only = TRUE)