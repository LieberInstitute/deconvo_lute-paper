#!/usr/bin/env R

#
# methods for calculating z
#
#

zsource_transform <- function(zsource, transformv){
  # transform zsource dataset
  #
  # zsource: matrix of data used to calculate z (e.g. counts matrix)
  #
}

zsource_type <- function(zsource, varname = "cellType_broad_hc"){
  # get the types (e.g. cell types.columns in z)
  #
  # zsource: matrix of data used to calculate z (e.g. counts matrix)
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

zsource_markerdata <- function(zsource, varname = "cellType_broad_hc", 
                           method = "mean_ratio", mr.assay = "logcounts"){
  # get the gene markers (e.g. features/rows in z)
  #
  # note:
  #   wraps several marker gene functions
  #
  # zsource: matrix of data used to calculate z (e.g. counts matrix)
  # varname : name of variable in zsource containing types
  # mr.assay : assay_name argument for get_mean_ratio2()
  #
  marker.data <- NA
  if(is(zsource, "SingleCellExperiment")){zsource <- counts(zsource)}
  message("checking varname exists in zsource")
  typev <- zsource_type(zsource = zsource, varname = varname)
  if(method == "mean_ratio"){
    if(!is(zsource, "SingleCellExperiment")){
      stop("method mean_ratio requires that zsource is a SingleCellExperiment")}
    message("using method mean_ratio")
    require(DeconvoBuddies)
    marker.data <- DeconvoBuddies::get_mean_ratio2(zsource, 
                                                   cellType_col = varname, 
                                                   assay_name = mr.assay,
                                                   add_symbol = TRUE)
  } else if(method == "findMarkers"){
    message("using method findMarkers")
    require(scran)
    # findMarkers code goes here
  }
  return(marker.data)
}

