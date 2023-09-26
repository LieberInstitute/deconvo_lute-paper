#!/usr/bin/env R

#
# Defines the cell size scale factors. Uses the following sets:
#
# * manul : the manually-defined set, also used in pseudobulk experiments.
# 
# * null : the null set, where factors are equal across types (both as "1").
#

# get sample s data from rnascope
rnascope.sizes <- metadata(sce.iter)[["cell.sizes"]]

list.sizes <- lapply(unique(rnascope.sizes$sample.id), function(sample.id){
  sizes.iter <- rnascope.sizes[rnascope.sizes$sample.id==sample.id,1]
  names(sizes.iter) <- rnascope.sizes[rnascope.sizes$sample.id==sample.id,2]
  sizes.iter <- sizes.iter[order(names(sizes.iter))]
  return(sizes.iter)
})

list.null <- lapply(unique(rnascope.sizes$sample.id), function(sample.id){
  c("glial" = 1, "neuron" = 1)
})

names(list.sizes) <- names(list.null) <- unique(rnascope.sizes$sample.id)

# get cell size factor series
#list.s.pred <- list(s.set.manual = c("glial" = 3, "neuron" = 10),
#                    s.set.null = c("glial" = 1, "neuron" = 1))
list.s.pred <- list(rnascope = list.sizes, null = list.null)