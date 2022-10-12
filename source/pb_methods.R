#!/usr/bin/env R

# 
# Methods for conducting experiments with random pseudobulking with some 
# cell-level reference data.
# 
# example:
# new_pseudobulk_expt <- get_pb_experiment()
#

#---------------------------------------
# pseudobulk experiment helper functions
#---------------------------------------
get_lpb <- function(sef, datv = NA, nj = NA, ctvarname = "celltype.treg", 
                    counts.summary.method = "mean", seed.num = 2,
                    scale.range = 500:2000, get.results = TRUE, lz = NA){
  # get list of pseudobulked counts tables
  #
  # arguments
  # sef : SummarizedExperiment, SingleCellExperiment, or similar, ideally 
  #   filtered on some marker genes or z features. Note, cell type labels in
  #   sef[[ctvarname]] should not contain "_" (underscore) characters.
  # datv : vector of relative cell type representation weights. Length should be 
  #   equal to nj*nk. Values here correspond to the mpb design matrix. If NA, this is 
  #   randomized and variable nj is used. For example, datv = c(1,2) implies the
  #   relative representation of c(0.33, 0.66) for the first and second types.
  # nj : number of pseudobulked samples to simulate. This is only used if datv==NA.
  # ctvarname : name of the cell type labels, contained in sef coldata.
  # counts.summary.method : method for summarizing counts to get y.data table.
  # seed.num : token for random seed, should be an integer.
  # scale.range: integer vector from which to randomly take total counts in
  #   sample-wise expression.
  # get.results : whether to make table of results.
  # 
  # returns
  # list of pseudobulk results, inc. option for results df.
  #
  # examples
  #
  require(dplyr)
  lpb <- list() # make return object
  if(is(sef, "SingleCellExperiment")){
    require(SingleCellExperiment);require(DelayedArray)
  }
  if(!ctvarname %in% colnames(colData(sef))){
    stop("variable ", ctvarname, " not in sef coldata.")}
  set.seed(seed.num); lpb <- list()
  # parse types
  kvarv <- sef[[ctvarname]]; klabv <- unique(kvarv); nk <- length(klabv)
  # parse pb data -- get datv using either arg datv,nj
  if(is(datv, "logical")){
    if(is(nj, "logical")){
      stop("provide either datv or nj.")
    } else{
      datv <- rep(sample(1e3, nk), nj)
    }
  }
  # get pseudobulk design matrix
  ncoln <- length(datv)/nk
  if(!ncoln %% 1 == 0){stop("invalid datv dimensions")}
  mpb <- matrix(datv, ncol = ncoln) %>% 
    apply(2, function(ci){ci/sum(ci)})
  rownames(mpb) <- klabv
  colnames(mpb) <- paste0("j_", seq(ncol(mpb)))
  lpb[["pi_pb"]] <- mpb
  # get sample scale factors, or randomized total counts
  scalev <- sample(scale.range, ncol(mpb))
  lpb[["scalev"]] <- scalev
  names(lpb[["scalev"]]) <- colnames(mpb)
  # set up counts for sampling
  ct <- assays(sef)$counts 
  ctlabv <- colnames(ct) <- paste0(kvarv, "_", seq(ncol(ct)))
  # get list of pseudobulked counts tables
  lct <- lapply(seq(ncol(mpb)), function(ji){
    # get sample-specific info
    scalej <- scalev[ji] # sample scale factor
    cellv.ij <- round(mpb[,ji]*scalej) # vector of total cell counts
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
  names(lct) <- colnames(mpb)
  lpb[["listed_counts_pb"]] <- lct
  if(!is.na(counts.summary.method)){
    if(counts.summary.method == "mean"){
      ypb <- do.call(cbind, lapply(lct, 
                                   function(ii){rowMeans(ii)}))
    } else if(counts.summary.method == "median"){
      ypb <- do.call(cbind, lapply(lct, 
                                   function(ii){rowMedians(ii)}))
    } else{
      stop("counts.summary.method not recognized")
    }
    colnames(ypb) <- colnames(mpb)
    rownames(ypb) <- rownames(sef)
    lpb[["y_data_pb"]] <- ypb
  }
  if(get.results & !is(lz, "logical")){
    df.res <- pb_report(lz.compare = lz, lpb = lpb)
    lpb[["pb_report"]] <- df.res
  }
  return(lpb)
}

get_pi_est <- function(z.data, y.data, method = "nnls", return.prop = TRUE){
  # get pi est using nnls
  #
  # method: type of strict deconvolution method to use (either "nnls", "glm" or 
  #   "bvls"). Defaults to "nnls".
  # return.prop: whether to return proportions. if False, returns magnitudes 
  #   from NNLS.
  #
  # example:
  # z.data <- matrix(sample(1000, 100), ncol = 5)
  # y.data <- matrix(sample(1000, 40), ncol = 2)
  # pi.est <- get_pi_est(z.data, y.data)
  # colSums(pi.est) # check within-sample types add to 1
  # # [1] 1 1
  #
  methodv.valid <- c("nnls", "glm", "bvls")
  if(method %in% methodv.valid){
    message("running ", method, "...")
    if(method == "nnls"){
      require(nnls)
      pi.dati <- do.call(rbind, lapply(seq(ncol(z.data)), 
                                       function(i){
                                         nnls::nnls(y.data, 
                                                    z.data[,i])$x}))
    } else if(method == "glm"){
      require(glmnet)
      pi.dati <- do.call(rbind, lapply(seq(ncol(z.data)), 
                                       function(i){
                                         coef(glmnet(y.data, 
                                                     z.data[,i], 
                                                     lambda = 0, 
                                                     lower.limits = 0, 
                                                     intercept = FALSE))}))
    } else{
      require(bvls)
      pi.dati <- do.call(rbind, lapply(seq(ncol(z.data)), 
                                       function(i){
                                         bvls(y.data, z.data[,ii],
                                              bl = rep(0, ncol(y.data)), 
                                              bu = rep(Inf, ncol(y.data)))$x
                                       }))
    }
  } else{
    stop("provide a valid method.")
  }
  if(return.prop){
    pi.dati <- apply(pi.dati, 2, function(ci){ci/sum(ci)}) 
  }
  colnames(pi.dati) <- colnames(y.data)
  rownames(pi.dati) <- colnames(z.data)
  return(pi.dati)
}

#----------------------------
# report data.frame functions
#----------------------------
# results df
pb_report <- function(lz.compare, lpb, method.str = "nnls", save.results = FALSE, 
                      save.fpath = "df-results_s-transform-expt_dlpfc-ro1.rda"){
  # makes results table for pseudobulk experiment.
  # 
  # note: 
  # 
  # * makes tall table of results comparing pi_est to pi_pb 
  #     (truth), similar to:
  # 
  # method1 cell_type sample_id pi_true pi_est diff_pb_minus_est
  # method2 cell_type sample_id pi_true pi_est diff_pb_minus_est
  #
  # arguments:
  # lz.compare : list of z tables to compare
  # cell.typev : ordered vector of cell type names
  # pi.pb : matrix of true cell proportions (cells X samples).
  # znamev : names of z tables to use in lz list.
  # save.results : whether to save final results table.
  # save.fpath : path to save final results table.
  #
  # example:
  # df.report <- pb_report(lz.compare=get_exe_lz(),lpb=get_exe_lpb())
  #
  require(reshape2)
  lz <- lz.compare
  znamev <- names(lz)
  if(!"pi_pb" %in% names(lpb)|
     !"y_data_pb" %in% names(lpb)|
     !"scalev" %in% names(lpb)){
    stop("lpb requires items 'pi_pb', 'y_data_pb', and 'scalev'.")
  }
  pi.pb.matrix <- 
    pi.pb <- 
    lpb[["pi_pb"]] 
  y.data <- lpb[["y_data_pb"]]
  cell.typev <- rownames(pi.pb)
  scalev <- lpb[["scalev"]]
  if(!ncol(lz[[1]])==length(cell.typev)){
    stop("num. k not equal in lz and lpb.")
  }
  df.tall <- do.call(rbind, 
                     lapply(seq(length(znamev)), 
                            function(ii){
                              znamei <- znamev[ii]
                              
                              message(znamei); z.data <- lz[[znamei]]
                              pi.dati <- get_pi_est(z.data, 
                                                    y.data, 
                                                    method = method.str)
                              # define cell type variables
                              pi.pb.df <- as.data.frame(pi.pb.matrix)
                              pi.dati <- as.data.frame(pi.dati)
                              pi.dati$cell_type <- 
                                pi.pb.df$cell_type <- 
                                cell.typev
                              # get tall tables
                              pi.dati.tall <- reshape2::melt(pi.dati, id = "cell_type")
                              pi.pb.tall <- reshape2::melt(pi.pb.df, id = "cell_type")
                              # return tall table
                              df.tall <- data.frame(cell_type = pi.dati.tall$cell_type,
                                                    sample_id = pi.dati.tall$variable,
                                                    pi_est = pi.dati.tall$value,
                                                    pi_true = pi.pb.tall$value,
                                                    pi_diff = pi.pb.tall$value-
                                                      pi.dati.tall$value)
                              df.tall$method <- znamei; return(df.tall)
                            }))
  # append scale info
  df.tall$scale <- NA
  for(ji in seq(length(unique(df.tall$sample_id)))){
    samplei <- unique(df.tall$sample_id)[ji]
    df.tall[df.tall$sample_id==samplei,]$scale <- scalev[ji]}
  if(save.results){save(df.tall, file = save.fpath)}
  return(df.tall)
}

#---------------
# plot functions
#---------------
# pi plot functions
pi_plot <- function(est, true){
  # pi_plot
  # 
  # est: vector of estimated pi values
  # true: vector of true pi values
  #
  #
  #
  require(ggplot2)
  est <- rnorm(1000, 1000, 500)/1000
  true <- rnorm(1000, 1000, 500)/1000
  ggpt <- ggplot(dfp, aes(x = true, y = est)) + theme_bw() +
    geom_point(alpha = 0.3) + xlab("pi_true") + ylab("pi_est") +
    geom_abline(intercept = 0, slope = 1, color = "red", 
                lwd = 1.2, alpha = 0.8) +
    geom_smooth(method = "lm", color = "blue",
                lwd = 1.2, alpha = 0.5)
  return(ggpt)
}

pi_plot_series <- function(df.tall = NA, ctvarname = NA, alpha.value = 0.4,
                           save.plots = TRUE, save.dpath = file.path(),
                           fn.handle = "newplots"){
  # makes plot series
  #
  # df.tall : valid tall report data.frame, e.g. such as returned by running 
  #   `pb_report()`.
  #
  # example
  # pi_plot_series()
  #
  # returns
  # lgg, list of ggplot objects
  #
  require(ggplot2)
  if(is(df.tall, "logical")){df.tall<-get_exe_dftall()}
  # parse variable names
  # get main plot object
  ggpt.main <- ggplot(df.tall, 
                      aes_string(x = pi_true, 
                          y = pi_est),
                      environment = environment()) +
    geom_abline(intercept = 0, 
                slope = 1, 
                color = "red", 
                lwd = 1.2, 
                alpha = alpha.value) +
    geom_smooth(method = "lm", 
                color = "blue", 
                lwd = 1.2, 
                alpha = alpha.value)
  # label series
  ggpt.all <- ggpt.main + 
    geom_point(alpha = alpha.value)
  ggpt.all.method <- ggpt.main + 
    geom_point(aes(color = method), 
               alpha = alpha.value)
  ggpt.all.sampleid <- ggpt.main + 
    geom_point(aes(color = sample_id), 
               alpha = alpha.value)
  ggpt.all.col <- ggpt.main + 
    geom_point(aes(color = cell_type), 
               alpha = alpha.value)
  # get facet series
  ggpt.all.celltype.facet <- ggpt.all.col + 
    facet_wrap(vars(cell_type), nrow = 1)
  ggpt.all.method.facet <- ggpt.all.col + 
    facet_wrap(~method, nrow = 1)
  ggpt.all.sampleid.facet <- ggpt.all.col + 
    facet_wrap(~sample_id, nrow = 1)
  # parse save options
  if(save.plots){
    message("saving new plots...")
    ggsave(paste0("ggpt-all_",fn.handle,".png"), 
           plot = ggpt.all, device = "png",
           path = save.dpath, width = 5, height = 5, units = "in", dpi = 400)
    ggsave(paste0("ggpt-all-method_",fn.handle,".png"), 
           plot = ggpt.all.method, device = "png",
           path = save.dpath, width = 5, height = 5, units = "in", dpi = 400)
    ggsave(paste0("ggpt-all-sampleid_",fn.handle,".png"), 
           plot = ggpt.all.sampleid, device = "png",
           path = save.dpath, width = 5, height = 5, units = "in", dpi = 400)
    ggsave(paste0("ggpt-all-col_",fn.handle,".png"), plot = ggpt.all.col, 
           device = "png",
           path = save.dpath, width = 5, height = 5, units = "in", dpi = 400)
    ggsave(paste0("ggpt-all-celltype-facet_",fn.handle,".png"), 
           plot = ggpt.all.celltype.facet, device = "png",
           path = save.dpath, width = 5, height = 5, units = "in", dpi = 400)
    ggsave(paste0("ggpt-all-method-facet_",fn.handle,".png"), 
           plot = ggpt.all.method.facet, device = "png",
           path = save.dpath, width = 5, height = 5, units = "in", dpi = 400)
    ggsave(paste0("ggpt-all-sampleid-facet_",fn.handle,".png"), 
           plot = ggpt.all.sampleid.facet, device = "png",
           path = save.dpath, width = 5, height = 5, units = "in", dpi = 400)
  }
  return(list(ggpt.main = ggpt.main, 
              ggpt.all = ggpt.all, 
              ggpt.all.method = ggpt.all.method,
              ggpt.all.sampleid = ggpt.all.sampleid,
              ggpt.all.col = ggpt.all.col, 
              ggpt.all.celltype.facet = ggpt.all.celltype.facet,
              ggpt.all.method.facet = ggpt.all.method.facet,
              ggpt.all.sampleid.facet = ggpt.all.sampleid.facet))
}

# scale plot functions
scale_plot_series <- function(df.tall = NA, alpha.value = 0.4, save.plots = TRUE,
                              save.dpath = file.path(), fn.handle = "newplots"){
  # makes plot series
  #
  #
  require(ggplot2)
  if(is(df.tall, "logical")){df.tall<-get_exe_dftall()}
  # pi est vs true, scale as size
  ggpt.scalesize <- ggplot(df.tall, 
                           aes(x = pi_true, 
                               y = pi_est, 
                               size = scale)) +
    geom_point(alpha = alpha.value)
  # get main template plot for series
  ggpt.main <- ggplot(df.tall, aes(x = scale))
  # get plot series
  ggpt.main.piest <- ggpt.main + 
    geom_point(aes(y = pi_est), 
               alpha = alpha.value)
  ggpt.main.pitrue <- ggpt.main + 
    geom_point(aes(y = pi_true), 
               alpha = alpha.value)
  ggpt.main.pidiff <- ggpt.main + 
    geom_point(aes(y = pi_diff), 
               alpha = alpha.value)
  ggpt.main.pidiff.sampleid.col <- ggpt.main + 
    geom_point(aes(y = pi_diff, 
                   color = sample_id), 
               alpha = alpha.value)
  # get facets
  ggpt.main.pidiff.sampleid.facet <- ggpt.main.pidiff + 
    facet_wrap(~sample_id)
  ggpt.main.pidiff.sampleid.celltype.facet <- ggpt.main.pidiff + 
    facet_wrap(~sample_id+cell_type)
  if(save.plots){
    message("saving new plots...")
    ggsave(paste0("ggpt-scalesize_",fn.handle,".png"), 
           plot = ggpt.scalesize, device = "png",
           path = save.dpath, width = 5, height = 5, units = "in", dpi = 400)
    ggsave(paste0("ggpt-scale-main-piest_",fn.handle,".png"), 
           plot = ggpt.main.piest, device = "png",
           path = save.dpath, width = 5, height = 5, units = "in", dpi = 400)
    ggsave(paste0("ggpt-scale-main-pitrue_",fn.handle,".png"), 
           plot = ggpt.main.pitrue, device = "png",
           path = save.dpath, width = 5, height = 5, units = "in", dpi = 400)
    ggsave(paste0("ggpt-scale-main-pidiff_",fn.handle,".png"), 
           plot = ggpt.main.pidiff, 
           device = "png",
           path = save.dpath, width = 5, height = 5, units = "in", dpi = 400)
    ggsave(paste0("ggpt-scale-main-pidiff-sampleid-col_",fn.handle,".png"), 
           plot = ggpt.main.pidiff.sampleid.col, device = "png",
           path = save.dpath, width = 5, height = 5, units = "in", dpi = 400)
  }
  return(list(
    ggpt.scalesize = ggpt.scalesize,
    ggpt.main = ggpt.main,
    ggpt.main.piest = ggpt.main.piest,
    ggpt.main.pitrue = ggpt.main.pitrue,
    ggpt.main.pidiff = ggpt.main.pidiff,
    ggpt.main.pidiff.sampleid.col = ggpt.main.pidiff.sampleid.col,
    ggpt.main.pidiff.sampleid.facet = ggpt.main.pidiff.sampleid.facet,
    ggpt.main.pidiff.sampleid.celltype.facet =
      ggpt.main.pidiff.sampleid.celltype.facet
  ))
}

#----------------
# example objects
#----------------
get_exe_lpb <- function(seed.num = 2){
  # get example lpb object
  #
  # example
  # get_exe_lpb()
  #
  require(SummarizedExperiment)
  set.seed(seed.num)
  # make counts data
  ct <- matrix(
    sample(100, 50*100, replace = T), 
    nrow = 50)
  # get summarized experiment
  sef <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = ct))
  sef[["celltypes"]] <- c(rep("excit", 40), 
                          rep("inhib", 30), 
                          rep("oligo", 10), 
                          rep("other", 20))
  # make pb series
  samp1.ratios <- c(10,10,5,10)
  samp2.ratios <- c(5,5,10,5)
  lpb <- get_lpb(sef, 
                 datv = c(samp1.ratios, 
                          samp2.ratios), 
                 ctvarname = "celltypes")
  return(lpb)
}

get_exe_lz <- function(seed.num = 2){
  # get example lz object
  #
  set.seed(seed.num)
  cell.typev = c("excit", "inhib", "oligo", "other")
  z.data <- matrix(sample(1000, 200), ncol = 4)
  colnames(z.data) <- paste0("k_", seq(ncol(z.data)))
  lz <- list(z1 = z.data, z2 = z.data)
  return(lz)
}

get_exe_sef <- function(ct.str = "celltype", seed.num = 2){
  # get example filtered summarized experiment
  #
  # ct.str : base string for cell types. should not include "_" (underscore).
  #
  #
  ct <- matrix(
    sample(100, 50*100, replace = T), 
    nrow = 50)
  sef <- SummarizedExperiment(assays = list(counts = ct))
  sef[["celltypes"]] <- paste0(rep("celltype", 100), seq(4))
  return(sef)
}

get_exe_dftall <- function(seed.num = 2){
  # get example tall report data
  #
  message("getting example df.tall data...")
  set.seed(seed.num)
  df.tall <- data.frame(cell_type = rep(c("k_1", "k_2"), 4),
                        sample_id = rep(paste0("j_", seq(2)), each = 2),
                        pi_est = rnorm(4, 100, 20)/100,
                        pi_true = rnorm(4, 100, 20)/100,
                        method = rep("z1", 4),
                        scale = c(rep(c(100,1000), each = 2)))
  df.tall$pi_diff <- df.tall$pi_true-df.tall$pi_est
  df.tall <- df.tall[,c(1:4, 7, 5:6)]
  return(df.tall)
}

#-----------------------------------
# main pseudobulk experiment wrapper
#-----------------------------------
get_pb_experiment <- function(lz = NA, scef = NA, 
                              datv = c(1,1,1,1), nj = NA,
                              scale.range = 500:2000,
                              ctvarname = "celltypes",
                              plot.results = TRUE,
                              method.str = "nnls", 
                              plot.save.dpath = "",
                              plot.fname.handle = "newplots",
                              save.plots = TRUE,
                              seed.num = 2){
  # run a pseudobulk experiment
  #
  # lz : list of z signature matrices for which to compare pi estimates.
  # scef : filtered SummarizedExperiment or similar object containing the gene
  #   markers or cell type-specific features.
  # datv : vector of relative cell type representation weights. Length should be 
  #   equal to nj*nk. Values here correspond to the mpb design matrix. If NA, this is 
  #   randomized and variable nj is used. For example, datv = c(1,2) implies the
  #   relative representation of c(0.33, 0.66) for the first and second types.
  # scale.range: range of total counts for random scaling of total expression.
  #   This is taken randomly for each simulated sample.
  # ctvarname : cell type variable name, should be contained in scef coldata.
  # plot.results : whether to make summary plots of various pi variables, scale
  #   versus pi differences, etc.
  # method.str : character string of strict deconvolution method to make pi_est,
  #   can be either "nnls" (default), "glm", or "bvls".
  # seed.num : integer for the random seed.
  # plot.save.dpath : path location to save new plots.
  # plot.fname.handle : character string to append to new plot filenames.
  # 
  # returns
  # randomiezed pseudobulk data, results table, optional plots
  #
  # example
  # pb.expt <- get_pb_experiment() # get experiment data
  # # plot results
  # attach(pb.expt$lpb$pb_report)
  # pb.expt$lgg_plots$lgg.scale$ggpt.scalesize
  # 
  lr <- list()
  if(is(sef, "logical")){sef <- get_exe_sef()}
  if(is(lz, "logical")){lz <- get_exe_lz()}
  lr[["lpb"]] <- get_lpb(sef = sef, lz = lz, datv = datv, nj = NA, 
                         ctvarname = ctvarname, seed.num = seed.num, 
                         scale.range = scale.range, get.results = TRUE)
  if(plot.results){
    lgg.pi <- pi_plot_series(lr[["pb_report"]], 
                             save.dpath = plot.save.dpath, 
                             fn.handle = plot.fname.handle,
                             save.plots = save.plots)
    lgg.scale <- scale_plot_series(lr[["pb_report"]],
                                   save.dpath = plot.save.dpath, 
                                   fn.handle = plot.fname.handle,
                                   save.plots = save.plots)
    lr[["lgg_plots"]] <- list("lgg.pi" = lgg.pi, "lgg.scale" = lgg.scale)
  }
  return(lr)
}
