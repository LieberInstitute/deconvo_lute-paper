#!/usr/bin/env R

#
#
#
#

libv <- c("SingleCellExperiment", "limma", "SummarizedExperiment", "scuttle")
sapply(libv, library, character.only = TRUE)

#---------------
# manage paths
#---------------
# get save dpath
code.dname <- "09_manuscript"
proj.dname <- "deconvo_method-paper"
save.dpath <- file.path(proj.dname, "outputs", code.dname)

#--------------
# manage params
#--------------
markeri <- "k2"
handle.str <- "ro1-dlpfc"

#----------
# load data
#----------
sce.fname <- paste0("sce_marker-adj-",markeri,"_",handle.str,".rda")
sce.fpath <- file.path(save.dpath, sce.fname)
scei <- get(load(sce.fpath))

#-----------------------------
# means and variances by donor
#-----------------------------
groupvar <- paste0(scei[["k2"]], ";", scei[["Sample"]])
table(groupvar)

scei[["groupvar"]] <- groupvar
sce.var <- aggregateAcrossCells(scei, ids = "groupvar", statistics = "var")

