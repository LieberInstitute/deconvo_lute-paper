#!/usr/bin/env R

# Author: Sean Maden
#
# Pre-filter and process bulk RNA-seq data.
#

# load data
source("deconvo_method-paper/code/11_bulk-expression-analysis/00_parameters.R")
sapply(libv, library, character.only = T)
rse <- get(load(file.path(rse.bulk.filepath)))
cd <- colData(rse)
# add batch id
cd[,batch.variable] <- paste0(cd[,donor.variable], "_", cd[,location.variable])
# add experiment groups
cond1 <- cd[,experiment.condition1]
cond2 <- cd[,experiment.condition2]
cd[,condition.variable] <- paste0(cond1, "_", cond2)
# reassign coldata
colData(rse) <- cd
# filter gene types
rd <- rowData(rse)
types.vector <- rd$gene_type
# save
save(rse, file = rse.bulk.path.new)
save(rse, file = rse.gene.filter.filepath)
