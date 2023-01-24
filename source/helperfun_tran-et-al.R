#!/usr/bin/env R
#
# Author: Sean Maden
#
# Helper functions for working with Tran et al 2021 multi-region brain snRNAseq
# datasets.
#
#

get_sce_tran <- function(loc.type = "dlpfc"){
  loc.str <- toupper(loc.type)
  addr <- paste0("https://libd-snrnaseq-pilot.s3.us-east-2.amazonaws.com/SCE_",
                 loc.str,"-n3_tran-etal.rda")
  fname <- unlist(strsplit(addr, "/"))
  fname <- fname[grepl("\\.rda$", fname)][1]
  if(!fname %in% list.files()){
    message("downloading SCE rda dataset...")
    utils::download.file(url = addr, destfile = fname)
  } else{message("found local SCE rda dataset")}
  message("loading SCE dataset...")
  return(get(load(fname)))
}

srcparam_tran <- function(){
  # source the key parameters for Tran et al 2021
  sce.celltype.varname <- "cellType"
  sce.donor.varname <- "donor"
  sce.k2.varname <- "k2"
  lparam <- list(varnames = list(celltype = sce.celltype.varname,
                                 donor = sce.donor.varname,
                                 klevels = list(k2 = sce.k2.varname)))
  return(lparam)
}

prepare_sce_tran <- function(expr.rescale.types = c("logcounts"), 
                             rescale.expr = FALSE, do.k.series = FALSE){
  #
  #
  # key steps to prepare downloaded sce object for analysis
  #
  
  lparam <- srcparam_tran()
  
  if(rescale.expr){
    if("logcounts" %in% expr.rescale.types){
      message("")
      if(!"logcounts" %in% names(assays(sce))){
        sce <- scuttle::logNormCounts(sce)
        save(sce, file = fname)
      }
    }
  }
  
  
  
  # make k variables from cell type variable
  if(do.k.series){
    kvarname <- "k2"
    if(!kvarname %in% colnames(colData(sce))){
      sce[[kvarname]] <- ifelse(grepl("Excit|Inhib", sce[[ctvarname]]), 
                                "neuron", "other")
      save(sce, file = fname)
    }
  }
  
  return(sce)
}