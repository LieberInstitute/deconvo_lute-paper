#!/usr/bin/env R

#
# Transformations for the z signature matrix.
#

s_rescale <- function(z, factorv = NULL, meanv = NULL, sdv = NULL,
                      disttype = c("norm"), seed.num = NULL, 
                      return.md = FALSE, constrain.nn = TRUE){
  # rescale z using some type-level factor s
  #
  # Notes: this type of transformation is used, for instance, to rescale cell
  #   types according to their expected sizes. 
  #
  # Note: if meanv and varv provided, select factor values randomly from
  # a distribution.
  # z : signature matrix (cols = k types, rows = g genes)
  # factorv : vector of rescale factors to be applied to types 
  #   (length == ncol(z))
  # meanv : vector of distribution means to be applied to types
  #   (length == ncol(z), order same as varv).
  # sd : vector of distribution variances, to be applied to types 
  #   (length == ncol(z), order same as meanv).
  # disttype : type of distribution, or distribution function, from which
  #   to randomly sample. Only applies if meanv and varv also provided.
  # seed.num : set the random seed. If NULL, is assigned using rnorm.
  # return.md : whether to return z_rescale with metadata.
  # constrain.nn : whether to constrain values in z_rescale to be non-negative.
  #   Changes negative values to 0 if TRUE.
  #
  # example:
  # s_rescale(z)
  z_rescale <- NA # newly rescaled z matrix output
  ms <- NA # matrix corresponding to the s factor values
  if(is.null(factorv)){
    if(is.null(meanv)|is.null(varv)){
      stop("need to provide factors or distribution info.")
    } else{
      if(length(meanv)==ncol(z) & length(sdv)==ncol(z)){
        if(disttype == "norm"){
          if(is.null(seed.num)){set.seed(seq(1000))}else{set.seed(seed.num)}
          # draw random transform values from dist
          ms <- do.call(rbind, 
                           lapply(seq(length(meanv)), function(ii){
            rnorm(n = nrow(z), mean = meanv[ii], sd = sdv[ii])
          }))
          z_rescale <- z*ms
        } else{
          stop("distribution function not recognized.")
        }
      } else{
        stop("lengths of meanv and varv must be equal to num. types in z.")
      }
    }
  } else{
    if(length(factorv)==ncol(z)){
      z_rescale <- z*rep(factorv, each = nrow(z))
      ms <- matrix(rep(factorv, each = nrow(z)), nrow = nrow(z))
    } else{
      stop("Num. factors in factorv must be equal to num. types in z.")}
  }
  if(constrain.nn){z_rescale[z_rescale<0] <- 0} # parse non-negativity
  lr <- z_rescale
  if(return.md){
    largs <- list(seed.num = seed.num, factorv = factorv, meanv = meanv, 
                  sdv = sdv, mfact = mfact, constrain.nn = constrain.nn)
    lmd <- list(s = ms, args = largs)
    lr <- list(z_rescale = z_rescale, metadata = lmd)
  }
  return(lr)
}



