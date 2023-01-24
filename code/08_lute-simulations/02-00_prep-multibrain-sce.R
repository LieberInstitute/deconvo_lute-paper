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