#!/usr/bin/env R

# Author: Sean Maden
#
# Get a series of quality control summary plots.
#

# load
source("deconvo_method-paper/code/11_bulk-expression-analysis/00_parameters.R")
sapply(libv, library, character.only = T)
rse.filter <- get(load(rse.gene.filter.filepath))

# set up data summaries
# rse data
cd <- colData(rse.filter)
assays.rse <- assays(rse.filter)
if(!assay.name %in% names(assays.rse)){
  stop("Error, assay.name ", assay.name, " not found in rse data.")} else{
    expression <- assays.rse[[assay.name]]
}

# get new cd variables
cd[,batch.variable] <- paste0(cd[,donor.id],"_",
                              cd[,brain.location.id])
cd[,condition.variable] <- paste0(cd[,experiment.condition1],"_",
                                  cd[,experiment.condition2])

# save new qc summary plots
for(type in type.vector){
  message("Working on summary type ", type, "...")
  get_summary_list(type = type, plot.fname = plot.fname, 
                   variable.vector = variable.vector, cd = cd, 
                   counts = expression, save.path = save.path)
  message("Finished with summary type ", type, ".")
}
message("Done with all summary  types.")
