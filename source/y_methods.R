#!/usr/bin/env R

# 
# Methods for the y mixed signals matrices. These have standard dimensions [G,J]
# for G gene marker features and J samples.
#

get_lpb <- function(scef, datv = NA, nj = NA, ctvarname = "celltype.treg", 
                    counts.summary.method = "mean", seed.num = 2,
                    scale.range = 500:2000, get.results = TRUE, lz = NA){
  # get list of pseudobulked counts tables
  #
  # arguments
  # datv : vector of relative cell type representation weights. Length should be 
  #   equal to nj*nk. Values here correspond to the mpb design matrix. If NA, this is 
  #   randomized and variable nj is used. For example, datv = c(1,2) implies the
  #   relative representation of c(0.33, 0.66) for the first and second types.
  # scef : SummarizedExperiment, SingleCellExperiment, or similar, ideally 
  #   filtered on some marker genes or z features.
  # nj : number of pseudobulked samples to simulate. This is only used if datv==NA.
  # counts.summary.method : method for summarizing counts to get y.data table.
  # scale.range: range of total counts for random scaling of total expression.
  #   This is taken randomly for each simulated sample.
  # get.results : whether to make table of results.
  # 
  # returns
  # list of pseudobulk results, inc. option for results df.
  #
  # examples
  #
  # require(SummarizedExperiment)
  # ct <- matrix(sample(100, 50*102, replace = T), nrow = 50)
  # sef <- SummarizedExperiment(assays = list(counts = ct))
  # sef[["celltypes"]] <- c(rep("glia", 2), rep("neuron", 100))
  # lpb <- get_lpb(sef, datv = c(1,10), ctvarname = "celltypes")
  # print(dim(lpb[[1]])) # [1] 2 1
  # print(dim(lpb[[2]][[1]])) # [1] 50 1474
  # print(dim(lpb[[3]])) # [1] 50  1
  # lpb[[1]]
  # #              j_1
  # #glia   0.09090909
  # #neuron 0.90909091
  # head(lpb[[3]])
  # #          j_1
  # #[1,] 54.36092
  # #[2,] 48.54138
  # #[3,] 42.43284
  # #[4,] 50.17843
  # #[5,] 52.94844
  # #[6,] 46.34600
  #
  require(dplyr)
  if(is(scef, "SingleCellExperiment")){
    require(SingleCellExperiment);require(DelayedArray)
  }
  set.seed(seed.num); lpb <- list()
  # parse types
  kvarv <- scef[[ctvarname]]; klabv <- unique(kvarv); nk <- length(klabv)
  # parse pb data -- get datv using either arg datv,nj
  if(is(datv, "logical")){
    if(is(nj, "logical")){
      stop("provide either datv or nj.")
    } else{
      datv <- rep(sample(1e3, nk), nj)
    }
  }
  ncol <- length(datv)/nk
  # get pseudobulk design matrix
  mpb <- matrix(datv, ncol=ncol) %>% apply(2, function(ci){ci/sum(ci)})
  rownames(mpb) <- klabv
  colnames(mpb) <- paste0("j_", seq(ncol(mpb)))
  lpb[["pi_pb"]] <- mpb
  scalev <- sample(scale.range, ncol(mpb)) # sample scale factors (total counts)
  ct <- assays(scef)$counts # set up counts for sampling
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
  names(lct) <- colnames(mpb)
  lpb[["scalev"]] <- scalev
  names(lpb[["scalev"]]) <- colnames(mpb)
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
    rownames(ypb) <- rownames(scef)
    lpb[["y_data_pb"]] <- ypb
  }
  if(get.results & !is(lz, "logical")){
    df.res <- pb_report(lz, cell.typev, pi.pb, 
                        znamev, save.results = F)
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

pb_report <- function(lz.compare, method.str = "nnls", 
                      save.results = FALSE, 
                      plot.results = TRUE,
                      cell.typev = c("Inhib", "Oligo", "other", "Excit"),
                      znamev = c("z.final", "zs1", "zs2"), 
                      pi.pb = lpb[["pi_pb"]],
                      save.fpath = "df-results_s-transform-expt_dlpfc-ro1.rda"){
  #
  #
  #
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
  #
  # example:
  #

require(SummarizedExperiment)

# get example data
method.str <- "nnls"
cell.typev = c("excit", "inhib", "oligo", "other")
z.data <- matrix(sample(1000, 200), ncol = 4)
colnames(z.data) <- paste0("k_", seq(ncol(z.data)))

# make counts data
ct <- matrix(sample(100, 50*100, replace = T), nrow = 50)

# get summarized experiment
sef <- SummarizedExperiment(assays = list(counts = ct))
sef[["celltypes"]] <- c(rep("excit", 40), rep("inhib", 30), 
                        rep("oligo", 10), rep("other", 20))
# make pb series
samp1.ratios <- c(10,10,5,10)
samp2.ratios <- c(5,5,10,5)
lpb <- get_lpb(sef, 
               datv = c(samp1.ratios, 
                        samp2.ratios), 
               ctvarname = "celltypes")

lz <- 
  lz.compare <- 
  list(z1 = z.data, 
       z2 = z.data)
pi.pb <- lpb$pi_pb
y.data <- lpb$y_data_pb

znamev <- names(lz)

require(reshape2)

lz <- lz.compare
znamev <- names(lz); 
pi.pb.matrix <- pi.pb 
y.data <- y.data
cell.typev <- rownames(pi.pb)

scalev <- lpb$scalev

df.tall <- do.call(rbind, lapply(seq(length(znamev)), function(ii){
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



  # 
  #
  require(reshape2)
  lz <- lz.compare; znamev <- names(lz); 
  pi.pb.matrix <- pi.pb 
  y.data <- y.data
  df.tall <- do.call(rbind, lapply(znamev, function(znamei){
    message(znamei); z.data <- lz[[znamei]]
    pi.dati <- get_pi_est(z.data, y.data, method = method.str)
    # get results table
    # define cell type variables
    pi.pb.df <- as.data.frame(pi.pb.matrix); pi.dati <- as.data.frame(pi.dati)
    pi.dati$cell_type <- pi.pb.df$cell_type <- cell.typev
    # get tall tables
    pi.dati.tall <- melt(pi.dati, id = "cell_type")
    pi.pb.tall <- melt(pi.pb.df, id = "cell_type")
    # return tall table
    df.tall <- data.frame(cell_type = pi.dati.tall$cell_type,
                          sample_id = pi.dati.tall$variable,
                          pi_est = pi.dati.tall$value,
                          pi_true = pi.pb.tall$value,
                          pi_diff = pi.pb.tall$value-pi.dati.tall$value)
    df.tall$method <- znamei; return(df.tall)
  }))
  if(save.results){save(df.tall, file = save.fpath)}
  return(df.tall)
}
