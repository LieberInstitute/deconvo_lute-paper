#!/usr/bin/env R

# Author: Sean Maden
#
# Pre-filter and process bulk RNA-seq data.
#

# load data
sapply(libv, library, character.only = T)
rse <- get(load(file.path(rse.bulk.filepath)))
# filter gene types
types.protein <- c("protein_coding")
types.nonpolya <- c("lncRNA", "Mt_rRNA", "rRNA", "Mt_tRNA")
types.include <- c(types.protein, types.nonpolya)
rd <- rowData(rse)
types.vector <- rd$gene_type
# save
rse.filename <- "rse_processed.rda"
save(rse, file = file.path(save.path, rse.filename))
