#!/usr/bin/env R

# Author: Sean Maden

#' musicParam-class
#'
#' Applies the MuSiC::music.basic() implementation of the MuSiC deconvolution 
#' algorithm.
#' 
#' @include lute_generics.R
#' @include referencebasedParam-class.R
#' 
#' @details Main constructor for class \linkS4class{musicParam}.
#' @rdname musicParam-class
#' @seealso 
#' \linkS4class{deconParam}
#' 
#' @examples
#' # example
#' lexample <- get_decon_example_data()
#' # param <- musicParam(s = lexample[["s"]], 
#' #                      y = lexample[["y"]], 
#' #                      z = lexample[["z"]])
#' 
#' # return only predicted proportions
#' # deconvolution(param)
#' 
#' # return full results
#' # param@return.info <- TRUE
#' # names(deconvolution(param))
#' 
#' @references 
#' 
#' Wang, Xuran and Jiaxin Fan. MuSiC: Multi-subject single cell deconvolution. 
#' (2022) GitHub, R package version 1.0.0. URL: https://github.com/xuranw/MuSiC.
#' 
#' Wang, X., Park, J., Susztak, K. et al. Bulk tissue cell type deconvolution 
#' with multi-subject single-cell expression reference. Nat Commun 10, 380 
#' (2019). https://doi.org/10.1038/s41467-018-08023-x
#' 
#' @aliases 
#' MuSiCParam-class
#'
setClass("musicParam", contains="referencebasedParam", 
         slots=c(sigma = "matrix", nu = "numeric", 
                 eps = "numeric", iter.max = "numeric"))

#' Make new object of class musicParam
#'
#' Main constructor for class \linkS4class{musicParam}.
#'
#' @param y Bulk mixed signals matrix of samples, which can be matched to single-cell samples.
#' @param z Signature matrix of cell type-specific signals. If not provided, can be computed from a
#' provided ExpressionSet containing single-cell data.
#' @param s Cell size factor transformations of length equal to the K cell types to deconvolve.
#' @param sigma Additional argument for algorithm.
#' @param nu Additional argument for algorithm.
#' @param iter.max Additional argument for algorithm.
#' @param eps Additional argument for algorithm.
#' @param return.info Whether to return metadata and original method outputs with predicted proportions.
#'
#' @examples
#' # example
#' lexample <- get_decon_example_data()
#' # param <- musicParam(s = lexample[["s"]], 
#' #                      y = lexample[["y"]], 
#' #                      z = lexample[["z"]])
#' 
#' # return only predicted proportions
#' # deconvolution(param)
#' 
#' # return full results
#' # param@return.info <- TRUE
#' # names(deconvolution(param))
#'
#' @returns New object of class \linkS4class{musicParam}.
#'
#' @details Takes standard inputs for the MuSiC algorithm
#' 
#' @export
musicParam <- function(y, z, s = NULL, sigma = NULL, nu = NULL, 
                       iter.max = NULL, eps = NULL, return.info = FALSE) {
  if(is(s, "NULL")){s <- rep(1, ncol(z))}
  if(is(sigma, "NULL")){sigma <- matrix(0, ncol = 1, nrow = nrow(z))}
  if(is(nu, "NULL")){nu <- 1e-10}
  if(is(iter.max, "NULL")){iter.max <- 1000}
  if(is(eps, "NULL")){eps <- 0}
  new("musicParam", y = y, z = z, s = s, sigma = sigma, nu = nu, 
      iter.max = iter.max, eps = eps, return.info = return.info)
}

#'
#'
#' @param list.pred Predictions list.
#' @param column.labels cell type labels.
#' @param row.labels sample id labels.
#'
#'
parse_deconvolution_predictions_results <- function(list.pred, 
                                                     column.labels, 
                                                     row.labels){
  if(is(column.labels, "NULL")){column.labels <- seq(length(list.pred[[1]]))}
  if(is(row.labels, "NULL")){row.labels <- seq(length(list.pred))}
  table.pred <- do.call(rbind, list.pred)
  table.pred <- apply(table.pred, 1, function(ri){ri/sum(ri)}) 
  table.pred <- t(table.pred)
  colnames(table.pred) <- column.labels
  rownames(table.pred) <- row.labels
  # convert
  table.pred <- cellProportionsPredictions(table.pred, column.labels, 
                                           row.labels)
  return(table.pred)
}

#' Deconvolution method for class \linkS4class{musicParam}
#' 
#' Main deconvolution method for the \linkS4class{musicParam} to run the 
#' \code{music.basic()} implementation of the MuSiC algorithm.
#' 
#' @param object An object of class \linkS4class{musicParam}.
#' 
#' @examples
#' # example
#' lexample <- get_decon_example_data()
#' # param <- musicParam(s = lexample[["s"]], 
#' #                      y = lexample[["y"]], 
#' #                      z = lexample[["z"]])
#' 
#' # return only predicted proportions
#' # deconvolution(param)
#' 
#' # return full results
#' # param@return.info <- TRUE
#' # names(deconvolution(param))
#' 
#' @returns Either a vector of predicted proportions, or a list containing 
#' predictions, metadata, and original outputs.
#'
#' @references 
#' 
#' Wang, Xuran and Jiaxin Fan. MuSiC: Multi-subject single cell deconvolution. 
#' (2022) GitHub, R package version 1.0.0. URL: https://github.com/xuranw/MuSiC.
#' 
#' Wang, X., Park, J., Susztak, K. et al. Bulk tissue cell type deconvolution 
#' with multi-subject single-cell expression reference. Nat Commun 10, 380 
#' (2019). https://doi.org/10.1038/s41467-018-08023-x
#'
#' @export
setMethod("deconvolution", signature(object = "musicParam"), function(object){
  lparam <- callNextMethod()
  # instantiate objects
  nu <- object[["nu"]]
  iter.max <- object[["iter.max"]]
  eps <- object[["eps"]]
  sigma <- object[["sigma"]]
  y <- lparam[["y"]]
  z <- lparam[["z"]]
  s <- lparam[["s"]]
  sigma <- as.matrix(sigma)
  bulk.samples.index.vector <- seq(ncol(y))
  result <- lapply(bulk.samples.index.vector, function(index){
    MuSiC::music.basic(X = z, 
                       Y = y[,index,drop=F], 
                       S = s, 
                       Sigma = sigma, 
                       nu = nu, 
                       iter.max = iter.max, 
                       eps = eps)
  })
  names(result) <- colnames(y)
  predictions <- lapply(result, function(iter){iter$p.weight})
  lr <- parse_deconvolution_predictions_results(predictions,
                                           colnames(z),
                                           names(predictions))
  if(object[["return.info"]]){
    lr <- list(predictions = predictions, 
               result.info = result, 
               metadata = lparam[["metadata"]])}
  return(lr)
})