#!/usr/bin/env R

#
# Get new MAE with nucleus area filtered RNAscope datasets.
#
#

repo.name <- "deconvo_method-paper"
folder.name <- "08_cell-counts-vs-qc-snrnaseq_notebook-scripts_dlpfc-cohort1"

#-----
# load
#-----
# load mae data
mae.filename <- "mae_with-rpkm_additional-data_final.rda"
mae <- get(load(file.path(repo.name, "outputs", "01_prepare-datasets", mae.filename)))

#-------
# source
#-------
# source mae utilties
script.path <- file.path(repo.name, "notebooks", folder.name, "mae_utilities.R")
source(script.path)

#-----
# main
#-----
# new data with filtered rnascope sce objects
mae <- new_mae_with_filters(mae)

# save
mae.filtered.name <- "mae_rnascope-filtered_cohort1.rda"
mae.filtered.path <- file.path(repo.name, "outputs", folder.name, mae.filtered.name)
save(mae, file = mae.filtered.path)
