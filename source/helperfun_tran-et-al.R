#!/usr/bin/env R
#
# Author: Sean Maden
#
# Helper functions for working with Multi-region Brain (a.k.a. mrb) datasets 
# published in Tran et al 2021.
#
#

get_sce_mrb <- function(check.fpath = ".", loc.type = "dlpfc", download = T,
                         http.addr = paste0("https://libd-snrnaseq-pilot",
                                            ".s3.us-east-2.amazonaws.com/")){
  # get_sce_tran
  #
  # check for sce, and optionally download if not found.
  #
  loc.str <- toupper(loc.type)
  sce.fname <- paste0("SCE_", loc.str,"-n3_tran-etal.rda")
  sce.fpath <- file.path(check.fpath, sce.fname)
  if(file.exists(sce.fpath)){
    message("found local SCE rda dataset")
    sce <- get(load(sce.fpath))
  } else{
    if(download){
      message("downloading sce...")
      addr <- paste0(http.addr, sce.fname)
      utils::download.file(url = addr, destfile = sce.fpath)
      sce <- get(load(sce.fpath))
    } else{
      message("sce not found at: ", sce.fpath)
      return(NULL)
    }
  }
  return(sce)
}

parse_k_brain <- function(sce, kvar = 2, celltypevar = "cellType"){
  cd <- colData(sce); ctvar <- cd[,celltypevar]
  kvar <- paste0("k", kvar)
  if(kvar %in% colnames(cd)){
    message("kvar ",kvar, " already in coldata.")
    return(sce)
  } else if(kvar == 2){
    newkvar <- ifelse(grepl("Excit|Inhib", ctvar), "neuron", "non-neuron")
    cd[,kvar] <- newkvar; colData(sce) <- cd
    return(sce)
  } else{
    stop("Error, invalid kvar.")
  }
  return(NULL)
}

prep_sce <- function(sce.fname = "SCE_DLPFC-n3_tran-etal.rda", 
                             rescale.assays = "logcounts", kvar = 2,
                             celltypevar = "cellType"){
  # prepare_sce_tran
  # 
  # key steps to prepare downloaded sce object for analysis
  #
  if("logcounts" %in% rescale.assays){
    if(!"logcounts" %in% names(assays(sce))){
      message("Getting logcounts...")
      sce <- scuttle::logNormCounts(sce)
    }
  }
  sce <- parse_k_brain(sce, kvar, celltypevar)
  return(sce)
}
