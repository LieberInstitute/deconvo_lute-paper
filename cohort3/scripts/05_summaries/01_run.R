#!/usr/bin/env R

# Author: Sean Maden
#
# Summarize dataset elements.
# 
# * "GSE107011_Processed_data_TPM.txt", GEO record.
#
#
#
#

libv <- c("SummarizedExperiment", "biomaRt")
sapply(libv, library, character.only = TRUE)
source("./source/00_read_experiment_data.R")
experimentData <- getExperimentData()
save.image("./env/05_summaries/01_run_script.RData")