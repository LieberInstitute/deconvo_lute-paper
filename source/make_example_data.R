#!/usr/bin/env R

#
# Makes some example data for use in a deconvolution task
#
# example:
# source("make_example_data.R")
# source("z_transform.R")
# ldecon <- ldecon_example()
#

zdata_example <- function(z.dist, k.value, z.nfeat, force.yz.nonneg, 
                          z.transformv, ltransform){
  # get example z data (e.g. feature values by type)
  #
  #
  #
  message("getting example z datasets...")
  z_original <- z_rescale <- NA
  # get z data
  z.dist.mean <- z.dist[1]; z.dist.var <- z.dist[2]
  z_original <- do.call(cbind, lapply(seq(k.value), function(i){
    z.dati <- rnorm(z.nfeat, z.dist.mean, z.dist.var)
    if(force.yz.nonneg){z.dati[z.dati<0] <- -1*z.dati[z.dati<0]}
    z.dati
  }))
  colnames(z_original) <- paste0("k_",seq(k.value))
  # do z transformations
  z.transformv <- z.transformv[z.transformv %in% c("s_rescale")]
  if(length(z.transformv) > 0){z_rescale <- z_original}
  if("s_rescale" %in% z.transformv){
    message("transforming z using `s_rescale` method...")
    lparam <- ltransform[["s_rescale"]]
    z_rescale <- s_rescale(z_rescale, factorv = lparam[["factorv"]],
                           constrain.nn = force.yz.nonneg)
  }
  if(is.na(z_rescale[1])){
    z_rescale <- z_original} # parse rescale default
  return(list(z_original = z_original, z_rescale = z_rescale))
}

ydata_example <- function(y.dist, y.nfeat, force.yz.nonneg, j.value = NA){
  # get example y data (e.g. feature values by sample)
  #
  #
  #
  message("getting example y datasets...")
  y_original <- y_rescale <- NA
  if(is.na(j.value)){
    j.value <- sample(100, 1)
  } else{if(j.value<1){stop("j.value must be > 0")}}
  y.dist.mean <- y.dist[1]; y.dist.var <- y.dist[2]
  y_original <- do.call(cbind, lapply(seq(j.value), function(i){
    y.dati <- rnorm(y.nfeat, y.dist.mean, y.dist.var)
    if(force.yz.nonneg){y.dati[y.dati<0] <- -1*y.dati[y.dati<0]}
    y.dati
  }))
  colnames(y_original) <- paste0("j_",seq(j.value))
  if(is.na(y_rescale)){y_rescale <- y_original} # parse rescale default
  return(list(y_original = y_original, y_rescale = y_rescale))
}

pi_example <- function(pi_data, pi_est, k.value, pi.est.funv, z.data, y.data, 
                       pi.method.validv = c("nnls")){
  # get example pi data (actual and estimated vectors of amounts by k type)
  #
  #
  #
  message("getting example pi values...")
  # get the pi amounts (proportions) by type
  if(is.na(pi_data)){pi_data <- rep(1/k.value, k.value)}
  # get pi estimates 
  # note: 
  #   this is the same as "strict deconvolution"
  #   filter on some list of valid methods, pi.method.validv
  # pi.method.validv <- c("nnls")
  if(is.na(pi_est)){skip.pi.est <- FALSE}
  if(skip.pi.est){
    message("using provided pi_est rather than calculating new pi_est.")
  } else{
    message("calculating pi_est using provided function(s)...")
    pi.est.funv <- tolower(pi.est.funv)
    pi.est.funv <- pi.est.funv[pi.est.funv %in% pi.method.validv]
    if(length(pi.est.funv)>0){
      pi_est <- lapply(pi.est.funv, function(funi){
        pi.dat <- NA
        if(funi=="nnls"){
          require(nnls)
          pi.dat <- do.call(rbind, lapply(seq(ncol(z.data)), function(i){
            nnls::nnls(y.data, z.data[,i])$x
          }))
          rownames(pi.dat) <- colnames(z.data)
          colnames(pi.dat) <- colnames(y.data)
        }
        pi.dat
      })
      names(pi_est) <- pi.est.funv
    }
  }
  return(list(pi_data = pi_data, pi_est = pi_est))
}

lmd_example <- function(z.dist.mean, z.dist.var, y.dist.mean, y.dist.var,
                        desc.str, seed.num, k.value, j.value, z.transformv,
                        source = "example"){
  # get metadata for example deconvolution datasets
  #
  #
  message("making example metadata...")
  z.data.info.misc <- list("dist.mean" = z.dist.mean, "dist.var" = z.dist.var)
  y.data.info.misc <- list("dist.mean" = y.dist.mean, "dist.var" = y.dist.var)
  y.data.info <- list(source = "example", misc = y.data.info.misc)
  z.data.info <- list(source = "example", misc = z.data.info.misc)
  desc.str <- "example random dataset"
  lmd <- list("description" = desc.str, "seed" = seed.num, 
              "k_value" = k.value, "j_value" = j.value,
              "z_info" = z.data.info, "y_info" = y.data.info,
              "z_transformations" = z.transformv)
  return(lmd)
}

ldecon_example <- function(seed.num = 0, k.value = 2, j.value = 10, 
                           n.feat = 1000, pi.data = NA, z.dist = c(5, 3),
                           y.dist = c(10, 5), force.yz.nonneg = TRUE,
                           z.transformv = c("s_rescale"),
                           ltransform = list(s_rescale = list(factorv = seq(2))),
                           pi.est.funv = c("nnls")){
  # example:
  # 
  # ldecon <- ldecon_example()
  #
  # seed.num : set the random seed. If NULL, is assigned using rnorm.
  # k.value : number of types/labels (e.g. cell types).
  # j.value: set the j total samples (e.g. num. bulk samples)
  # n.feat : number of features for z and y (e.g. num. genes)
  # pi.data : vector of type proportions. used to make pseudobulk. should have
  #   length == k.value or number of types. if null, vector of equal proportions
  #   is generated.
  # z.dist: vector for z normal distribution (mean, var)
  # y.dist: vector for y normal distribution (mean, var)
  # force.yz.nonneg : whether to force y and z data to be non-negative
  # z.transformv : vector of transformations to perform for z.
  # ltransform : list of transformation function arguments. Names in list correspond
  #   to valid transformation functions.
  # pi.est.funv : vector of functions to use to estimate the pi vector of type
  #   amounts (e.g. cell type proportions). If "nnls", uses `nnls::nnls()`.
  #
  # returns:
  # ldecon, list of deconvolution experiment entities
  #
  #
  # set value defaults
  # z_original <- z_rescale <- NA; y_original <- y_rescale <- NA
  pi_data <- pi_est <- lmd <- NA
  set.seed(seed.num) # get seed info
  # get the k total types 
  # (e.g. cell types)
  if(is.na(k.value)){
    k.value <- sample(c(2:10), 1)
  } else{if(k.value < 2){stop("k.value must be > 1")}}
  # get pi, z, y
  z.nfeat <- y.nfeat <- n.feat # features the same in z and y
  # get z data
  lz <- zdata_example(z.dist = z.dist, k.value = k.value, z.nfeat = z.nfeat, 
                      force.yz.nonneg = force.yz.nonneg, 
                      z.transformv = z.transformv, 
                      ltransform = ltransform)
  # get y data
  ly <- ydata_example(y.dist = y.dist, y.nfeat = y.nfeat, 
                      force.yz.nonneg = force.yz.nonneg, 
                      j.value = j.value)
  # get pi data
  lpi <- pi_example(pi_data = pi_data, pi_est = pi_est, k.value = k.value, 
                    pi.est.funv = pi.est.funv, z.data = lz[[2]], y.data = ly[[2]], 
                    pi.method.validv = c("nnls"))
  # parse metadata list
  lmd <- lmd_example(z.dist.mean = z.dist[1], z.dist.var = z.dist[2], 
                     y.dist.mean = y.dist[1], y.dist.var = y.dist[2],
                     desc.str = desc.str, seed.num = seed.num, 
                     k.value = k.value, j.value = j.value, 
                     z.transformv = z.transformv,
                     source = source)
  # return ldecon
  ldecon <- list("pi_data" = lpi[["pi_data"]], "pi_est" = lpi[["pi_est"]], 
                 "y_original" = ly[["y_original"]], 
                 "y_rescale" = ly[["y_rescale"]],
                 "z_original" = lz[["z_original"]], 
                 "z_rescale" = lz[["z_rescale"]],
                 "metadata" = lmd)
  return(ldecon)
}