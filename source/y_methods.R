#!/usr/bin/env R

# 
# Methods for the y mixed signals matrices. These have standard dimensions [G,J]
# for G gene marker features and J samples.
#

get_lpb <- function(scef, datv = NA, nj = NA, ctvarname = "celltype.treg", 
                    counts.summary.method = "mean", seed.num = 2){
  # get list of pseudobulked counts tables
  #
  # datv : vector of cell type scale data. Length should be equal to nj*nk. 
  #   Values here correspond to the mpb design matrix. If NA, this is randomized 
  #   and variable nj is used.
  # scef : SingleCellExperiment, ideally filtered on some marker genes or z features.
  # nj : number of pseudobulked samples to simulate. This is only used if datv==NA.
  # counts.summary.method : method for summarizing counts to get y.data table.
  # 
  set.seed(seed.num)
  # parse types
  kvarv <- scef[[ctvarname]]; klabv <- unique(kvarv); nk <- length(klabv)
  # parse pb data -- get datv using either arg datv,nj
  if(is.na(datv)){
    if(is.na(nj)){
      stop("provide either datv or nj.")
    } else{
      datv <- rep(sample(1e3, nk), nj)
    }
  }
  ncol <- length(datv)/nk
  # get pseudobulk design matrix
  mpb <- matrix(datv, ncol=ncol) %>% apply(2, function(ci){ci/sum(ci)})
  rownames(mpb) <- klabv; colnames(mpb) <- paste0("j_", seq(ncol(mpb)))
  scalev <- sample(1000:10000, ncol(mpb)) # sample scale factors (total counts)
  ct <- counts(scef) # set up counts for sampling
  ctlabv <- colnames(ct) <- paste0(kvarv, "_", seq(ncol(ct)))
  # get list of pseudobulked counts tables
  lct <- lapply(seq(ncol(mpb)), function(ji){
    # get sample-specific info
    scalej <- scalev[ji] # sample scale factor
    cellv.ij <- mpb[,ji]*scalej # vector of total cell counts
    # get randomized counts data
    ct.pb.j <- do.call(cbind, lapply(klabv, function(ki){
      num.cells.ij <- cellv.ij[ki] # num cells to sample for this type
      cnvf <- ctlabv[grepl(ki, gsub("_.*", "", ctlabv))]
      cnvf.index.ij <- cnvf[sample(seq(length(cnvf)), 
                                   num.cells.ij, replace = T)]
      return(ct[,cnvf.index.ij])
    }))
    return(ct.pb.j)
  })
  if(counts.summary.method == "mean"){
    ypb <- do.call(cbind, lapply(lct, function(ii){rowMeans(ii)}))
  } else if(counts.summary.method == "median"){
    ypb <- do.call(cbind, lapply(lct, function(ii){rowMedians(ii)}))
  } else{
    stop("counts.summary.method not recognized.")
  }
  rownames(ypb) <- names(lct)
  lpb <- list("listed_counts_pb" = lct, "y_data_pb" = ypb, "pi_pb" = mpb)
  return(lpb)
}
