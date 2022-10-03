#!/usr/bin/env R

#
# methods for calculating z
#
#

zsource_transform <- function(zsource, transformv){
  # transform zsource dataset
  #
  #
  #
}

zsource_type <- function(zsource, varname = "cellType_broad_hc"){
  # get the types (e.g. cell types.columns in z)
  #
  #
  #
  typev <- NA
  if(is(zsource, "SingleCellExperiment")){
    typev <- zsource[[varname]]} else{
      if(varname %in% colnames(zsource)){zsource[,varname]} else{
        stop("varname not a column name in zsource.")
      }
    }
  return(typev)
}

zsource_marker <- function(zsource, markerv, method = "mean_ratio"){
  # get the gene markers (e.g. features/rows in z)
  #
  #
  #
  markerv <- NA
  if(method == "mean_ratio"){require(DeconvoBuddies)}
  if(is(zsource, "SingleCellExperiment")){zsource <- counts(zsource)}
  
}

