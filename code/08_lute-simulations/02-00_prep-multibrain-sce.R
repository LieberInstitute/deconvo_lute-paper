#!/usr/bin/env R

# Author: Sean Maden 
#
# Prepping multiregion brain dataset for analysis.

libv <- c("SingleCellExperiment", "SummarizedExperiment")

# manage paths
code.dname <- "08_lute-simulations"
proj.dname <- "deconvo_method-paper"
save.dpath <- file.path(proj.dname, "outputs", code.dname)
source.dpath <- file.path(proj.dname, "source")
script.fname <- "helperfun_tran-et-al.R"
source(file.path(source.dpath, script.fname))

# get data
sce.fname <- "sce-mrb_dlpfc.rda"
sce.fpath <- file.path(save.dpath, sce.fname)
sce <- get_sce_mrb()
sce <- prep_sce()
save(sce, file = sce.fpath)

# check attributes
names(assays(sce))
# [1] "counts"    "logcounts"
table(sce[["k2"]])
# neuron  other 
# 3968   7234

# get markers
mr <- get_mean_ratio2(sce, "k2", "logcounts")
mr.fname <- "markers-k2_db-mr2_sce-dlpfc-mrb.rda"
mr.fpath <- file.path(save.dpath, mr.fname)
save(mr, file = mr.fpath)
