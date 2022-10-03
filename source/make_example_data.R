#!/usr/bin/env R

#
# Makes some example data for use in a deconvolution task
#

ldecon_example <- function(seed.num = NA, k.value = 2,
                           j.value = 10, n.feat = 1000,
                           pi.data = NA, z.dist = c(5, 3),
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
  # get seed info
  if(is.na(seed.num)){seed.num <- rnorm(1,100,100)}
  set.seed(seed.num)
  # get the k total types 
  # (e.g. cell types)
  if(is.na(k.value)){
    k.value <- sample(c(2:10), 1)
  } else{if(k.value < 2){stop("k.value must be > 1")}}
  # 
  # get pi, z, y
  # features the same in z and y
  z.nfeat <- y.nfeat <- n.feat 
  # get the pi amounts (proportions) by type
  if(is.na(pi.data)){pi.data <- rep(1/k.value, k.value)}
  # get z data
  z.dist.mean <- z.dist[1]; z.dist.var <- z.dist[2]
  z.data <- do.call(cbind, lapply(seq(k.value), function(i){
    z.dati <- rnorm(z.nfeat, z.dist.mean, z.dist.var)
    if(force.yz.nonneg){z.dati[z.dati<0] <- -1*z.dati[z.dati<0]}
    z.dati
  }))
  colnames(z.data) <- paste0("k_",seq(k.value))
  # do z transformations
  if("s_rescale" %in% z.transformv){
    message("transforming z using `s_rescale` method...")
    lparam <- ltransform[["s_rescale"]]
    z.data <- s_rescale(z.data, factorv = lparam[["factorv"]],
                        constrain.nn = force.yz.nonneg)
  }
  # get y data
  if(is.na(j.value)){
    j.value <- sample(100, 1)
  } else{if(j.value<1){stop("j.value must be > 0")}}
  y.dist.mean <- y.dist[1]; y.dist.var <- y.dist[2]
  y.data <- do.call(cbind, lapply(seq(j.value), function(i){
    y.dati <- rnorm(y.nfeat, y.dist.mean, y.dist.var)
    if(force.yz.nonneg){y.dati[y.dati<0] <- -1*y.dati[y.dati<0]}
    y.dati
  }))
  colnames(y.data) <- paste0("j_",seq(j.value))
  # get pi estimates 
  # note: this is the same as "strict deconvolution"
  pi.est <- NA
  # note: filter on some list of valid methods, pi.method.validv
  pi.method.validv <- c("nnls")
  pi.est.funv <- tolower(pi.est.funv)
  pi.est.funv <- pi.est.funv[pi.est.funv %in% pi.method.validv]
  if(length(pi.est.funv)>0){
    pi.est <- lapply(pi.est.funv, function(funi){
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
    names(pi.est) <- pi.est.funv
  }
  # parse metadata list
  z.data.info.misc <- list("dist.mean" = z.dist.mean,
                           "dist.var" = z.dist.var)
  y.data.info.misc <- list("dist.mean" = y.dist.mean,
                           "dist.var" = y.dist.var)
  y.data.info <- list(source = "example",
                      misc = y.data.info.misc)
  z.data.info <- list(source = "example",
                           misc = z.data.info.misc)
  desc.str <- "example random dataset"
  lmd <- list("description" = desc.str, "seed" = seed.num, 
              "k_value" = k.value, "j_value" = j.value,
              "z_info" = z.data.info, "y_info" = y.data.info,
              "z_transformations" = z.transformv)
  # return ldecon
  ldecon <- list("pi_data" = pi.data, "y_data" = y.data, 
                 "z_data" = z.data, "pi_est" = pi.est, 
                 "metadata" = lmd)
  return(ldecon)
}