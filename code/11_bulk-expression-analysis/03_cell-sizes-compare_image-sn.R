#!/usr/bin/env R

# Author: Sean Maden
#
# Compare cell sizes from snRNA-seq and image analysis.

source("deconvo_method-paper/code/07_cell-size-estimates/00_parameters.R")
sapply(libv, library, character.only = T)
sn.total <- get(load(sce.sizes.total.expression.path))
sn.gene <- get(load(sce.sizes.expressed.genes.path))
image.all <- get(load(image.cell.sizes.save.path))

# format cell size data
lsize <- list(total.expression = sn.total$cell.size.summary,
              expressed.genes = sn.gene$cell.size.summary,
              image.area = image.all$Nucleus_Area,
              marker.expression = image.all$AKT3_Copies)
dfsize <- do.call(rbind, lapply(seq(length(lsize)), function(index){
  size.type.label <- names(lsize)[index]
  dfsize.index <- lsize[[index]]
  dfsize.return <- data.frame(cell.type = dfsize.index[,1],
                              value = dfsize.index[,2])
  dfsize.return$size.type <- size.type.label
  return(dfsize.return)
}))
# get k2 labels
dfsize$k2 <- ifelse(dfsize[,1] %in% terms.neuron, "neuron",
                    ifelse(dfsize[,1] %in% terms.glial, "glial", "other"))
# get k2 sizes table
by.list <- list(k2 = dfsize[,"k2"], size.type = dfsize[,"size.type"])
dfsize.k2 <- aggregate(dfsize[,"value"], by = by.list, FUN = "median")
colnames(dfsize.k2)[3] <- "value"

# save
save(dfsize, file = table.cell.size.path)
save(dfsize.k2, file = table.cell.size.k2.path)
