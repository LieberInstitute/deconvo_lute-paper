#!/usr/bin/env R

# Author: Sean Maden
#
# Pre-filter and process bulk RNA-seq data.
#

# load data
source("deconvo_method-paper/code/11_bulk-expression-analysis/00_parameters.R")
sapply(libv, library, character.only = T)
rse <- get(load(file.path(rse.bulk.filepath)))
# add experiment groups
cd <- colData(rse)
cond1 <- cd[,experiment.condition1]
cond2 <- cd[,experiment.condition2]
colData(rse)[,condition.variable] <- paste0(cond1, "_", cond2)
# save rse to outputs
save(rse, file = rse.bulk.path.new)
# filter gene types
rd <- rowData(rse)
types.vector <- rd$gene_type
# save
save(rse, file = file.path(save.path, rse.gene.filter.filename))
