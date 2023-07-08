#!/usr/bin/env R

# Author: Sean Maden
#
# Pre-filter and process bulk RNA-seq data.
#

# load data
source("deconvo_method-paper/code/11_bulk-expression-analysis/00_parameters.R")
sapply(libv, library, character.only = T)
rse <- get(load(file.path(rse.bulk.filepath)))
# view gene types
rd <- rowData(rse)
table(rd$gene_type)
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
gene.types.vector <- rd$gene_type
gene.type.filter <- which(gene.types.vector %in% gene.types.include) 
rse.filtered <- rse[gene.type.filter,]
# save
save(rse, file = rse.bulk.path.new)
save(rse.filtered, file = rse.gene.filter.filepath)
