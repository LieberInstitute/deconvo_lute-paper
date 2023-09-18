#!/usr/bin/env R

# Author: Sean Maden
#
# Get differentially expressed genes (DEGs) among bulk sample experiment conditions.
#

source("deconvo_method-paper/code/11_bulk-expression-analysis/00_parameters.R")
sapply(libv, library, character.only = T)
rse.markers <- get(load(rse.k2markers.filepath))

# compare mean expression by gene across groups
expression.data <- assays(rse.markers)[[assay.name]] %>% 
  as.data.frame() %>% t() %>% as.data.frame()
expression.data$group <- colData(rse.markers)[,condition.variable]
expression.summary <- aggregate(expression.data, 
                                by = list(group = expression.data$group),
                                FUN = "mean")
expression.summary <- expression.summary[,c(1:ncol(expression.summary)-1)]

# compare mean expression differences by gene, across groups
groups <- expression.summary[,1]
genes.index <- seq(2,ncol(expression.summary))
expression.differences <- do.call(cbind, lapply(groups, function(group){
  groups.compare <- groups[!groups==group]
  filter <- expression.summary[,1]==group
  expression.group <- expression.summary[filter, genes.index] %>% as.matrix()
  expression.compare <- expression.summary[!filter, genes.index] %>% as.matrix()
  expression.difference <- sweep(expression.compare, 1, expression.group, FUN = "-")
  rownames(expression.difference) <- paste0(group, "-", groups.compare)
  expression.difference
}))

