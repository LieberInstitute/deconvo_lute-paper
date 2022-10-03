#!/usr/bin/env R

#
# methods for calculating z from zsource, where z is the features-by-types matrix
# e.g. genes-by-cell types or [G,K], and zsource is the data from which z is
# calculated (e.g. a counts matrix or a SingleCellExperiment from which counts 
# are obtainable).
#
#

zsource_transform <- function(zsource, transformv){
  # transform zsource dataset
  #
  # zsource: matrix of data used to calculate z (e.g. counts matrix)
  #
}

zsource_type <- function(zsource, type.varname = "cellType_broad_hc"){
  # get the types (e.g. cell types.columns in z)
  #
  # zsource: matrix of data used to calculate z (e.g. counts matrix)
  # type.varname: name of the types variable, identifiable from zsource.
  #
  typev <- NA
  if(is(zsource, "SingleCellExperiment")){
    typev <- zsource[[type.varname]]} else{
      if(type.varname %in% colnames(zsource)){
        zsource[,type.varname]} else{
        stop("type.varname not a column name in zsource.")
      }
    }
  return(typev)
}

zsource_markerdata <- function(zsource, type.varname = "cellType_broad_hc", 
                           method = "mean_ratio", mr.assay = "logcounts"){
  # get the gene markers (e.g. features/rows in z)
  #
  # note:
  #   wraps several marker gene functions
  #
  # zsource: matrix of data used to calculate z (e.g. counts matrix)
  # varname : name of variable in zsource containing types
  # method : name of method to identify marker genes from zsource.
  # mr.assay : assay_name argument for get_mean_ratio2()
  #
  # returns
  #
  # marker.data : data containing the features of interest and the values by 
  #   type (e.g. cell type) in zsource.
  #
  markerdata <- NA
  message("checking type.varname exists in zsource")
  typev <- zsource_type(zsource = zsource, 
                        type.varname = type.varname)
  if(method == "mean_ratio"){
    if(!is(zsource, "SingleCellExperiment")){
      stop("method mean_ratio requires that zsource is a SingleCellExperiment")}
    message("using method mean_ratio")
    require(DeconvoBuddies)
    markerdata <- DeconvoBuddies::get_mean_ratio2(zsource, 
                                                  cellType_col = type.varname, 
                                                  assay_name = mr.assay,
                                                  add_symbol = TRUE)
  } else if(method == "findMarkers"){
    message("using method findMarkers")
    if(is(zsource, "SingleCellExperiment")){zsource <- counts(zsource)}
    require(scran)
    # findMarkers code goes here
  }
  return(markerdata)
}

get_z_experiment <- function(ngenes.byk = NA, type.summary = "mean",
                             zsource = NA, type.varname = "cellType_broad_hc", 
                             method = "mean_ratio", mr.assay = "logcounts", 
                             markerdata = NA, type = NA, transform = NA){
  # calculate z for a deconvolution experiment
  #
  # zsource: matrix of data used to calculate z (e.g. counts matrix)
  # type.varname : name of variable in zsource containing types
  # method : name of method to identify marker genes from zsource.
  # mr.assay : assay_name argument for get_mean_ratio2()
  # marker.data : data containing the features of interest and the values by 
  #   type (e.g. cell type) in zsource.
  # ngenes.byk: top genes to select by type. If NA, retains all marker genes 
  #   provided upstream.
  # type.summary: function to use to summarize data by type.
  # 
  require(dplyr)
  if(is.na(markerdata)){
    markerdata <- zsource_markerdata(zsource, type.varname, method, mr.assay)
  }
  # get types
  typev <- unique(zsource[[type.varname]])
  # get z
  # get marker data for top genes
  ma <- markerdata
  top.markers <- lapply(typev, function(ki){
    markerv <- NA
    if(method == "mean_ratio"){
      # note: parses output from DeconvoBuddies::get_mean_ratio2()
      colnames(ma)[2] <- "celltype.target"
      mi <- ma %>% filter(celltype.target == ki) %>% arrange(rank_ratio)
      if(ngenes.byk < nrow(mi)){markerv <- markerv[seq(ngenes.byk)]}
    }
    message("keeping ",length(markerv)," features for type ",ki)
    markerv})
    
    
    
    # get k summaries
    if(type.summary == "mean"){
      ki.summary <- rowMeans(markerdati[,ki.summary.cols])
      rownames(ki.summary) <- ki.markers
    }
    ki.summary
  })))
  # summarize top genes
  md.byk <- do.call(cbind, lapply(seq(length(markerdata)), function(ii)))
}

