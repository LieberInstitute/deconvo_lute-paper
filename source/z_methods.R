#!/usr/bin/env R

#
# methods for calculating z from zsource, where z is the features-by-types matrix
# e.g. genes-by-cell types or [G,K], and zsource is the data from which z is
# calculated (e.g. a counts matrix or a SingleCellExperiment from which counts 
# are obtainable).
#
#
# example:
# sce <- get(load(sce.fpath))
# z.expt <- get_z_experiment(sce)
#

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

get_top_markers <- function(markerdata, typev, ngenes.byk, 
                            method.markers = "mean_ratio"){
  # get top markers from markerdata
  #
  # markerdata:
  # typev:
  # ngenes.byk: top markers to return for each cell type
  # method.markers: method used to generate markerdata.
  #
  # returns:
  # top markers as vector
  #
  ma <- markerdata
  ma.top <- do.call(rbind, lapply(typev, function(ki){
    mi <- NA
    if(method.markers == "mean_ratio"){
      # note: parses output from DeconvoBuddies::get_mean_ratio2()
      colnames(ma)[2] <- "celltype.target"
      mi <- ma %>% filter(celltype.target == ki) %>% 
        arrange(rank_ratio)
      if(nrow(mi)>0 & ngenes.byk<nrow(mi)){
        mi <- mi[seq(ngenes.byk),]}
    }
    message("keeping ",nrow(mi)," features for type ",ki); mi}))
  return(ma.top)
}

get_z_experiment <- function(zsource, 
                             markerdata = NA,
                             method.markers = "mean_ratio", 
                             mr.assay = "logcounts",
                             ngenes.byk = 20, 
                             type.varname = "cellType_broad_hc", 
                             summary.varname = "donor",
                             k.summary.method = "mean",
                             z.summary.method = "mean",
                             save.dpath = "deconvo_method-paper/outputs/",
                             marker.plots = TRUE,
                             return.all = FALSE){
  # calculate z for a deconvolution experiment
  #
  # zsource: data (of type SingleCellExperiment) used to calculate z (e.g. counts matrix)
  # markerdata : data containing the features of interest and the values by 
  #   type (e.g. cell type) in zsource.
  # method.markers : name of method to identify marker genes from zsource.
  # mr.assay : assay_name argument for get_mean_ratio2()
  # type.varname : name of variable in zsource coldata, containing type labels
  # summary.varname : name of variable in zsource coldata, containing labels on
  #   which to summarize the zsource counts.
  # k.summary.method : name of method for summarizing zsource counts on variable
  #   summary.varname and k types.
  # z.summary.method : name of method for summarizing within the k types for z
  # save.dpath : path of dir location to save new files (e.g. rds, png, etc.)
  # marker.plots : whether to generate new marker plots.
  # return.all: whether to return all intermediate outputs. If FALSE, returns 
  #   only the final z table.
  #
  # return:
  #   data for the z tables, with option to include intermediate outputs (see 
  #   argument return.all)
  #
  # example:
  # get_z_experiment(sce)
  # 
  require(dplyr)
  #if(!is(zsource, "SingleCellExperiment")){
  #  stop("zsource should be of type SingleCellExperiment")}
  if(is(markerdata, "logical")){
    markerdata <- zsource_markerdata(zsource = zsource, 
                                     type.varname = type.varname, 
                                     method = method.markers, 
                                     mr.assay = mr.assay)
  }
  # get variable levels
  typev <- unique(zsource[[type.varname]])
  summaryv <- unique(zsource[[summary.varname]])
  message("filtering source data on top markers")
  ma.top <- get_top_markers(ma = markerdata, typev = typev, 
                            ngenes.byk = ngenes.byk, 
                            method = method.markers)
  zsf <- zsource[rownames(zsource) %in% ma.top$gene,]
  # summarize counts by type and by category (e.g. donors)
  message("summarizing counts by type and ",summary.varname, "...")
  zs <- do.call(cbind, lapply(summaryv, function(ii){
    si <- do.call(cbind, lapply(typev, function(jj){
      zsf.filt <- zsf[[type.varname]] == jj & 
        zsf[[summary.varname]] == ii
      dati <- NA
      if(k.summary.method == "mean"){
        dati <- rowMeans(counts(zsf[,zsf.filt]))}
      dati
    }))
    colnames(si) <- paste0(ii, ";", typev);si
  }))
  # get z by summarizing type by category (e.g. donors)
  message("summarizing counts by type for final z table..."); cnv.zs <- colnames(zs)
  z <- t(apply(zs, 1, function(ri){
    unlist(lapply(typev, function(ki){
      dati <- NA; datf <- ri[grepl(ki, names(ri))]
      if(length(datf) > 0){
        if(z.summary.method == "mean"){
          dati <- mean(datf, na.rm = T)
        }
      }
      dati
    }))
  }))
  colnames(z) <- typev; lr <- z
  if(return.all){
    lr <- list(top.marker.data = ma.top, z.summary.filt = zs, z.final = z)
    if(marker.plots){
      # source("z_figures.R")
      lr[["plots"]] <- get_lgg_markers(df.markers = ma.top, save.dpath = save.dpath)
    }
  }
  return(lr)
}

