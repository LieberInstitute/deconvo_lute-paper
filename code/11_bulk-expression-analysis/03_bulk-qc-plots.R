#!/usr/bin/env R

# Author: Sean Maden
#
# Get a series of quality control summary plots.
#

# load
source("deconvo_method-paper/code/11_bulk-expression-analysis/00_parameters.R")
sapply(libv, library, character.only = T)
rse.filter <- get(load(rse.gene.filter.filepath))
cd <- colData(rse.filter)
assays.rse <- assays(rse.filter)
if(!assay.name %in% names(assays.rse)){
  stop("Error, assay.name ", assay.name, " not found in rse data.")} else{
    expression <- assays.rse[[assay.name]]
}
# save new qc summary plots
for(type in type.vector){
  message("Working on summary type ", type, "...")
  plot.filename <- paste0("ggjitter-composite_qc-rse_",
                          gsub("\\.", "-", type),".jpg")
  get_summary_list(type = type, plot.filename = plot.filename, 
                   variable.vector = variable.vector, 
                   cd = cd, expression = expression, 
                   save.path = save.path)
  message("Finished with summary type ", type, ".")
}
message("Done with all summary  types.")
