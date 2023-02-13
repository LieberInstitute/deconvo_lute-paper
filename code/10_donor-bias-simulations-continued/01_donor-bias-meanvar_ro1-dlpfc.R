#!/usr/bin/env R

#
#
#
#

libv <- c("SingleCellExperiment", "limma", "SummarizedExperiment")
sapply(libv, library, character.only = TRUE)

#---------------
# manage paths
#---------------
# get save dpath
code.dname <- "09_manuscript"
proj.dname <- "deconvo_method-paper"
save.dpath <- file.path(proj.dname, "outputs", code.dname)

# params
markeri <- "k2"
handle.str <- "ro1-dlpfc"

#----------
# load data
#----------
sce.fname <- paste0("sce_marker-adj-",markeri,"_",handle.str,".rda")
sce.fpath <- file.path(save.dpath, sce.fname)
scei <- get(load(sce.fpath))
