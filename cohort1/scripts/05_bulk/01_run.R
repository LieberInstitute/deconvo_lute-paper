#!/usr/bin/env R

# Author: Sean Maden
#
# Perform bulk conditions A/B tests.
#
#

libv <- c("scuttle")
sapply(libv, library, character.only = TRUE)

#---------
# load mae
#---------
new.mae.filename <- "mae_allsamples.rda"
mae.final.filepath <- file.path("outputs", "01_mae", new.mae.filename)
mae <- get(load(mae.final.filepath))

#-------------
# run a/b test
#-------------
base.path <- file.path("scripts/05_bulk/abtest/")
# prep dependencies
source(file.path(base.path, "00_helper-functions.R"))
# prep data
source(file.path(base.path, "01-01_prepare-mae.R")) # prep mae data
source(file.path(base.path, "01-02_prep-s-vector.R"))
# run experiments
source(file.path(base.path, "02-01-01_rse-counts_counts-yz_shared-reference-experiments.R"))
source(file.path(base.path, "02-01-02_rse-counts_logcounts-lutearg_shared-reference-experiments.R"))
source(file.path(base.path, "02-02-01_rse-counts_counts-yz_within-reference-experiments.R"))
source(file.path(base.path, "02-02-02_rse-counts_lognorm-yz_within-reference-experiments.R"))
source(file.path(base.path, "02-03-01_rse-rpkm_counts-yz_shared-reference-experiments.R"))
source(file.path(base.path, "02-03-02_rse-rpkm_logcounts-lutearg_shared-reference-experiments.R"))
source(file.path(base.path, "02-04-01_rse-rpkm_counts-yz_within-reference-experiments.R"))
source(file.path(base.path, "02-04-02_rse-rpkm_lognorm-yz_within-reference-experiments.R"))
source(file.path(base.path, "03_prep-experiment-results.R"))
###