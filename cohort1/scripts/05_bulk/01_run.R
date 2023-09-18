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

#------------------------
# filter training samples
#------------------------
# sample.id.train <- get(load("./outputs/00_preprocess/list_snrnaseq_sampleid.rda"))[["train"]]

# remove validation samples
dim(mae[["bulk.rnaseq"]])
cd.mae <- colData(mae)
cd.mae$sample.id.new <- gsub("_.*", "", cd.mae$sample.id)

validation.sample.id <- c("Br6432", "Br6522", "Br8667")
filter.string <- paste0(validation.sample.id, collapse = "|")
filter.mae <- !grepl(filter.string, cd.mae$sample.id.new)
table(filter.mae)
# filter.mae
# FALSE  TRUE 
# 7    15 

mae <- mae[,filter.mae,]
dim(mae[["bulk.rnaseq"]])

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

# save environment
save.image(file = "./env/05_bulk/01_run_script.RData")