#!/usr/bin/env R

# Author: Sean Maden
#
# Main parameters, or dependency objects, for deconvolution framework trials.

# dependencies
libv <- c("lute", "dplyr", "ggplot2", "gridExtra", "SingleCellExperiment", 
          "SummarizedExperiment", "scran")
sapply(libv, library, character.only = TRUE)

# save path
save.path <- file.path("deconvo_method-paper", "outputs", "02_pseudobulk-predictions")

#-------------------
# dlpfc cohort2 data
#-------------------
# mrb sce path
sce.mrb.path <- file.path("deconvo_method-paper", "outputs", "01_prepare-datasets", "sce-mrb_dlpfc.rda")

#-------------------
# dlpfc cohort1 data
#-------------------
# dlpfc markers path
sce.markers.list.path <- file.path("deconvo_method-paper", "outputs", "01_prepare-datasets", "list-scef_markers-k2-k3-k4_ro1-dlpfc.rda")

# define experiment function
order_svector <- function(s.vector = c("neuron" = 10, "glial" = 3)){
  s.vector[names(s.vector)[order(names(s.vector))]]
}

get_ypb_experiment_series_shuffle <- function(sce, sample.id.variable = "Sample", 
                                      celltype.variable = "k2", assay.name = "logcounts",
                                      sample.id.s.scale.shuffle.reference = NULL,
                                      s.data.frame = data.frame("glial" = 3, "neuron" = 10),
                                      algorithm.name = "nnls", return.dimensions = c("wide", "tall"),
                                      dfp.tall.errors = TRUE, system.sleep.sec = 2,
                                      deconvolution.algorithm = "nnls"){
  # s.vector <- order_svector(s.vector)
  # get pseudobulk experiment series, testing cellsize adjustment
  # use with dfp_tall_by_celltype()
  # get experiment results
  if(!sample.id.variable %in% colnames(colData(sce))){stop("Error: couldn't find sample.id.variable in sce coldata.")}
  celltype.variable.format <- sce[[celltype.variable]] %>% as.factor()
  unique.celltypes <- unique(celltype.variable.format)
  num.celltypes <- length(unique(celltype.variable.format))
  
  # get s vector pred for shuffle
  s.sample.id.vector <- unique(s.data.frame$sample.id)
  if(is(sample.id.s.scale.shuffle.reference, "NULL")){
    sample.id.pred.all <- s.sample.id.vector[1]
  } else{
    sample.id.pred.all <- sample.id.s.scale.shuffle.reference
  }
  message("using sample.id ",sample.id.pred.all, " for s shuffle reference.")
  s.vector.pred.all <- s.data.frame[s.data.frame$sample.id == sample.id.pred.all,
                                    grepl("s\\..*", colnames(s.data.frame))]
  s.vector.pseudobulk.pred.all <- as.numeric(s.vector.pred.all)
  names(s.vector.pseudobulk.pred.all) <- gsub("s\\.", "", names(s.vector.pred.all))
  #
  s.vector.null <- rep(1, num.celltypes); names(s.vector.null) <- unique.celltypes
  s.vector.null <- s.vector.null[order(match(names(s.vector.null), names(s.vector.pseudobulk.pred.all)))]
  
  # iterate on pseudobulk sample ids (don't change s.pred vector ....)
  s.sample.id.vector <- s.sample.id.vector[s.sample.id.vector %in% unique(sce[[sample.id.variable]])]
  dfp <- do.call(rbind, lapply(s.sample.id.vector, function(sample.id.iter){
    s.vector.iter <- s.data.frame[s.data.frame$sample.id == sample.id.iter,
                                  grepl("s\\..*", colnames(s.data.frame))]
    s.vector.pseudobulk.iter <- as.numeric(s.vector.iter)
    names(s.vector.pseudobulk.iter) <- gsub("s\\.", "", names(s.vector.iter))
    
    # get experiment proportions
    dfp.noscale <- get_ypb_experiment_results(sce[,sce[[sample.id.variable]]==sample.id.iter], 
                                              sample.id.variable = sample.id.variable, 
                                              celltype.variable = celltype.variable, 
                                              assay.name = assay.name,
                                              s.vector.ypb = s.vector.pseudobulk.iter,
                                              s.vector.pred = s.vector.null,
                                              system.sleep.sec = system.sleep.sec,
                                              deconvolution.algorithm = deconvolution.algorithm)
    
    dfp.withscale <- get_ypb_experiment_results(sce[,sce[[sample.id.variable]]==sample.id.iter], 
                                                sample.id.variable = sample.id.variable, 
                                                celltype.variable = celltype.variable, 
                                                assay.name = assay.name,
                                                s.vector.ypb = s.vector.pseudobulk.iter,
                                                s.vector.pred = s.vector.pseudobulk.pred.all,
                                                system.sleep.sec = system.sleep.sec,
                                                deconvolution.algorithm = deconvolution.algorithm)
    dfp.noscale$s.sample.id.pred.all <- sample.id.pred.all
    dfp.withscale$s.sample.id.pred.all <- sample.id.pred.all
    if(return.dimensions == "tall"){
      # get plot data -- tall
      dfp.noscale$type <- 'noscale'
      dfp.withscale$type <- 'withscale'
      dfp.noscale$sample.id <- rownames(dfp.noscale)
      dfp.withscale$sample.id <- rownames(dfp.withscale)
      dfp <- rbind(dfp.noscale, dfp.withscale)
      # if(dfp.tall.errors){dfp <- dfp_tall_append_errors(dfp)}
    } else{
      # get plot data -- wide
      colnames(dfp.noscale) <- paste0(colnames(dfp.noscale), ".noscale")
      colnames(dfp.withscale) <- paste0(colnames(dfp.withscale), ".withscale")
      identical(rownames(dfp.noscale), rownames(dfp.withscale))
      dfp <- cbind(dfp.noscale, dfp.withscale)
    }
    return(dfp)
  }))
  return(dfp)
}


get_ypb_experiment_series_shuffle_zpb <- function(sce.foriter, sce.deconvo, sample.id.variable = "Sample", 
                                      celltype.variable = "k2", assay.name = "logcounts",
                                      s.vector = c("glial" = 3, "neuron" = 10),
                                      algorithm.name = "nnls", return.dimensions = c("wide", "tall"),
                                      dfp.tall.errors = TRUE, system.sleep.sec = 2,
                                      deconvolution.algorithm = "nnls"){
  #
  # shuffles the z pseudobulk, using the same reference for deconvolution runs.
  # 
  # sce.foriter : sce containing sample groups for pb iterations.
  # sce.deconvo : the main sce object to use for all deconvolutions as the Z reference.
  #
  #

  s.vector <- order_svector(s.vector)
  # get pseudobulk experiment series, testing cellsize adjustment
  # use with dfp_tall_by_celltype()
  # get experiment results
  if(!sample.id.variable %in% colnames(colData(sce))){stop("Error: couldn't find sample.id.variable in sce coldata.")}
  celltype.variable.format <- sce[[celltype.variable]] %>% as.factor()
  unique.celltypes <- unique(celltype.variable.format)
  num.celltypes <- length(unique(celltype.variable.format))
  s.vector.null <- rep(1, num.celltypes); names(s.vector.null) <- unique.celltypes
  s.vector.null <- s.vector.null[order(match(names(s.vector.null), names(s.vector)))]
  
  # get experiment proportions
  # sce.iter <- sce[,sce[[sample.id.variable]]==sample.id.iter]
  dfp.noscale <- get_ypb_experiment_results_shuffle_zpb(sce.foriter, sce.deconvo, 
                                             sample.id.variable = sample.id.variable, 
                                             celltype.variable = celltype.variable, 
                                             assay.name = assay.name,
                                             s.vector.ypb = s.vector,
                                             s.vector.pred = s.vector.null,
                                          system.sleep.sec = system.sleep.sec,
                                          deconvolution.algorithm = deconvolution.algorithm)
  dfp.withscale <- get_ypb_experiment_results_shuffle_zpb(sce.foriter, sce.deconvo, 
                                               sample.id.variable = sample.id.variable, 
                                               celltype.variable = celltype.variable, 
                                               assay.name = assay.name,
                                               s.vector.ypb = s.vector,
                                               s.vector.pred = s.vector,
                                            system.sleep.sec = system.sleep.sec,
                                            deconvolution.algorithm = deconvolution.algorithm)
  if(return.dimensions == "tall"){
    # get plot data -- tall
    dfp.noscale$type <- 'noscale'
    dfp.withscale$type <- 'withscale'
    dfp.noscale$sample.id <- rownames(dfp.noscale)
    dfp.withscale$sample.id <- rownames(dfp.withscale)
    dfp <- rbind(dfp.noscale, dfp.withscale)
    if(dfp.tall.errors){dfp <- dfp_tall_append_errors(dfp)}
  } else{
    # get plot data -- wide
    colnames(dfp.noscale) <- paste0(colnames(dfp.noscale), ".noscale")
    colnames(dfp.withscale) <- paste0(colnames(dfp.withscale), ".withscale")
    identical(rownames(dfp.noscale), rownames(dfp.withscale))
    dfp <- cbind(dfp.noscale, dfp.withscale)
  }
  return(dfp)
}


get_ypb_experiment_series_shuffle_zs <- function(sce.foriter, sce.deconvo, sample.id.variable = "Sample", 
                                      celltype.variable = "k2", assay.name = "logcounts",
                                      sample.id.s.scale.shuffle.reference = NULL,
                                      s.data.frame = data.frame("glial" = 3, "neuron" = 10),
                                      algorithm.name = "nnls", return.dimensions = c("wide", "tall"),
                                      dfp.tall.errors = TRUE, system.sleep.sec = 2,
                                      deconvolution.algorithm = "nnls"){
  #
  # Shuffle the ZS reference across predictions
  #
  # sce.foriter : sce object to source pseudobulks. iterates on sample id groups when specified.
  # sce.deconvo : sce to use for all deconvolution experiments.
  #
  #
  #
  # s.vector <- order_svector(s.vector)
  # get pseudobulk experiment series, testing cellsize adjustment
  # use with dfp_tall_by_celltype()
  # get experiment results
  #
  #
  #
  if(!sample.id.variable %in% colnames(colData(sce))){stop("Error: couldn't find sample.id.variable in sce coldata.")}
  celltype.variable.format <- sce[[celltype.variable]] %>% as.factor()
  unique.celltypes <- unique(celltype.variable.format)
  num.celltypes <- length(unique(celltype.variable.format))
  
  # get s vector pred for shuffle
  s.sample.id.vector <- unique(s.data.frame$sample.id)
  if(is(sample.id.s.scale.shuffle.reference, "NULL")){
    sample.id.pred.all <- s.sample.id.vector[1]
  } else{
    sample.id.pred.all <- sample.id.s.scale.shuffle.reference
  }
  message("using sample.id ",sample.id.pred.all, " for s shuffle reference.")
  s.vector.pred.all <- s.data.frame[s.data.frame$sample.id == sample.id.pred.all,
                                    grepl("s\\..*", colnames(s.data.frame))]
  s.vector.pseudobulk.pred.all <- as.numeric(s.vector.pred.all)
  names(s.vector.pseudobulk.pred.all) <- gsub("s\\.", "", names(s.vector.pred.all))
  #
  s.vector.null <- rep(1, num.celltypes); names(s.vector.null) <- unique.celltypes
  s.vector.null <- s.vector.null[order(match(names(s.vector.null), names(s.vector.pseudobulk.pred.all)))]
  






  # iterate on pseudobulk sample ids (don't change s.pred vector ....)
  s.sample.id.vector <- s.sample.id.vector[s.sample.id.vector %in% unique(sce[[sample.id.variable]])]
  dfp <- do.call(rbind, lapply(s.sample.id.vector, function(sample.id.iter){
    s.vector.iter <- s.data.frame[s.data.frame$sample.id == sample.id.iter,
                                  grepl("s\\..*", colnames(s.data.frame))]
    s.vector.pseudobulk.iter <- as.numeric(s.vector.iter)
    names(s.vector.pseudobulk.iter) <- gsub("s\\.", "", names(s.vector.iter))
    
    # get experiment proportions
    sce.foriter <- sce[,sce[[sample.id.variable]]==sample.id.iter]
    dfp.noscale <- get_ypb_experiment_results_shuffle_z(sce.foriter = sce.foriter, 
                                              sce.deconvo = sce.deconvo,
                                              sample.id.variable = sample.id.variable, 
                                              celltype.variable = celltype.variable, 
                                              assay.name = assay.name,
                                              s.vector.ypb = s.vector.pseudobulk.iter,
                                              s.vector.pred = s.vector.null,
                                              system.sleep.sec = system.sleep.sec,
                                              deconvolution.algorithm = deconvolution.algorithm)
    
    dfp.withscale <- get_ypb_experiment_results_shuffle_z(sce.foriter = sce.foriter, 
                                                sce.deconvo = sce.deconvo,
                                                sample.id.variable = sample.id.variable, 
                                                celltype.variable = celltype.variable, 
                                                assay.name = assay.name,
                                                s.vector.ypb = s.vector.pseudobulk.iter,
                                                s.vector.pred = s.vector.pseudobulk.pred.all,
                                                system.sleep.sec = system.sleep.sec,
                                                deconvolution.algorithm = deconvolution.algorithm)
    dfp.noscale$s.sample.id.pred.all <- sample.id.pred.all
    dfp.withscale$s.sample.id.pred.all <- sample.id.pred.all
    if(return.dimensions == "tall"){
      # get plot data -- tall
      dfp.noscale$type <- 'noscale'
      dfp.withscale$type <- 'withscale'
      dfp.noscale$sample.id <- rownames(dfp.noscale)
      dfp.withscale$sample.id <- rownames(dfp.withscale)
      dfp <- rbind(dfp.noscale, dfp.withscale)
      # if(dfp.tall.errors){dfp <- dfp_tall_append_errors(dfp)}
    } else{
      # get plot data -- wide
      colnames(dfp.noscale) <- paste0(colnames(dfp.noscale), ".noscale")
      colnames(dfp.withscale) <- paste0(colnames(dfp.withscale), ".withscale")
      identical(rownames(dfp.noscale), rownames(dfp.withscale))
      dfp <- cbind(dfp.noscale, dfp.withscale)
    }
    return(dfp)
  }))
  return(dfp)
}

get_ypb_experiment_results <- function(sce, sample.id.variable = "Sample", 
                                       celltype.variable = "k2", assay.name = "logcounts",
                                       s.vector.ypb = c("glial" = 3, "neuron" = 10),
                                       s.vector.pred = c("glial" = 1, "neuron" = 1),
                                       deconvolution.algorithm = "nnls", system.sleep.sec = 2){
  s.vector.ypb <- order_svector(s.vector.ypb)
  s.vector.pred <- order_svector(s.vector.pred)
  # get results for a single iteration of an experiment
  # use with get_ypb_experiment_series()
  if(assay.name == "logcounts" & !"logcounts" %in% names(assays(sce))){sce <- scuttle::logNormCounts(sce)}
  unique.sample.id.vector <- unique(sce[[sample.id.variable]])
  dfp <- do.call(rbind, lapply(unique.sample.id.vector, function(sample.id){
    sce.iter <- sce[,sce[[sample.id.variable]]==sample.id]
    ypb.iter <- ypb_from_sce(sce = sce.iter, assay.name = assay.name, 
                             celltype.variable = celltype.variable,
                             sample.id.variable = sample.id.variable,
                             S = s.vector.ypb) %>% as.matrix()
    prop.true.iter <- table(sce.iter[[celltype.variable]]) %>% prop.table() %>% as.matrix() %>% t()
    prop.pred.iter <- lute(sce = sce.iter, y = ypb.iter, assay.name = assay.name, 
                           celltype.variable = celltype.variable, s = s.vector.pred, 
                           typemarker.algorithm = NULL, return.info = FALSE,
                           deconvolution.algorithm = deconvolution.algorithm)$deconvolution.results@predictions.table
    colnames(prop.pred.iter) <- paste0(colnames(prop.pred.iter), ".pred")
    colnames(prop.true.iter) <- paste0(colnames(prop.true.iter), ".true")
    dfp.iter <- cbind(prop.true.iter, prop.pred.iter) %>% as.data.frame()
    Sys.sleep(system.sleep.sec)
    dfp.iter
  }))
  rownames(dfp) <- unique.sample.id.vector
  return(dfp)
}

get_ypb_experiment_results_shuffle_zs <- function(sce.foriter, sce.deconvo, sample.id.variable = "Sample", 
                                       celltype.variable = "k2", assay.name = "logcounts",
                                       s.vector.ypb = c("glial" = 3, "neuron" = 10),
                                       s.vector.pred = c("glial" = 1, "neuron" = 1),
                                       deconvolution.algorithm = "nnls", system.sleep.sec = 2){
  #
  # use z shuffling
  #
  # sce.foriter : sce on which to base pseudobulks. iterates on sample id groups where specified.
  # sce.deconvo : sce to use across all deconvolution experiments (not used in pseudobulks directly)
  #
  s.vector.ypb <- order_svector(s.vector.ypb)
  s.vector.pred <- order_svector(s.vector.pred)
  # get results for a single iteration of an experiment
  # use with get_ypb_experiment_series()
  sce <- sce.foriter
  if(assay.name == "logcounts" & !"logcounts" %in% names(assays(sce))){sce <- scuttle::logNormCounts(sce)}
  unique.sample.id.vector <- unique(sce[[sample.id.variable]])
  dfp <- do.call(rbind, lapply(unique.sample.id.vector, function(sample.id){
    sce.iter <- sce[,sce[[sample.id.variable]]==sample.id]
    ypb.iter <- ypb_from_sce(sce = sce.iter, assay.name = assay.name, 
                             celltype.variable = celltype.variable,
                             sample.id.variable = sample.id.variable,
                             S = s.vector.ypb) %>% as.matrix()
    prop.true.iter <- table(sce.iter[[celltype.variable]]) %>% prop.table() %>% as.matrix() %>% t()
    prop.pred.iter <- lute(sce = sce.deconvo, y = ypb.iter, assay.name = assay.name, 
                           celltype.variable = celltype.variable, s = s.vector.pred, 
                           typemarker.algorithm = NULL, return.info = FALSE,
                           deconvolution.algorithm = deconvolution.algorithm)$deconvolution.results@predictions.table
    colnames(prop.pred.iter) <- paste0(colnames(prop.pred.iter), ".pred")
    colnames(prop.true.iter) <- paste0(colnames(prop.true.iter), ".true")
    dfp.iter <- cbind(prop.true.iter, prop.pred.iter) %>% as.data.frame()
    Sys.sleep(system.sleep.sec)
    dfp.iter
  }))
  rownames(dfp) <- unique.sample.id.vector
  return(dfp)
}

get_ypb_experiment_results_shuffle_zpb <- function(sce.foriter, sce.deconvo, sample.id.variable = "Sample", 
                                       celltype.variable = "k2", assay.name = "logcounts",
                                       s.vector.ypb = c("glial" = 3, "neuron" = 10),
                                       s.vector.pred = c("glial" = 1, "neuron" = 1),
                                       deconvolution.algorithm = "nnls", system.sleep.sec = 2){
  #
  # use z pseudobulk shuffling
  #
  # sce.foriter : sce on which to base pseudobulk iterations. iterates on sample id groups where specified.
  # sce.deconvo : sce to use across all deconvolution experiments (not used in pseudobulks directly)
  #

  s.vector.ypb <- order_svector(s.vector.ypb)
  s.vector.pred <- order_svector(s.vector.pred)
  # get results for a single iteration of an experiment
  # use with get_ypb_experiment_series()
  if(assay.name == "logcounts" & !"logcounts" %in% names(assays(sce.deconvo))){sce.deconvo <- scuttle::logNormCounts(sce.deconvo)}
  if(assay.name == "logcounts" & !"logcounts" %in% names(assays(sce.foriter))){sce.foriter <- scuttle::logNormCounts(sce.foriter)}

  unique.sample.id.vector <- unique(sce.foriter[[sample.id.variable]])
  dfp <- do.call(rbind, lapply(unique.sample.id.vector, function(sample.id){
    message(sample.id)
    sce.iter <- sce.foriter[,sce.foriter[[sample.id.variable]]==sample.id]
    ypb.iter <- ypb_from_sce(sce = sce.iter, assay.name = assay.name, 
                             celltype.variable = celltype.variable,
                             sample.id.variable = sample.id.variable,
                             S = s.vector.ypb) %>% as.matrix()
    prop.true.iter <- table(sce.iter[[celltype.variable]]) %>% prop.table() %>% as.matrix() %>% t()
    prop.pred.iter <- lute(sce = sce.deconvo, y = ypb.iter, assay.name = assay.name, 
                           celltype.variable = celltype.variable, s = s.vector.pred, 
                           typemarker.algorithm = NULL, return.info = FALSE,
                           deconvolution.algorithm = deconvolution.algorithm)$deconvolution.results@predictions.table
    colnames(prop.pred.iter) <- paste0(colnames(prop.pred.iter), ".pred")
    colnames(prop.true.iter) <- paste0(colnames(prop.true.iter), ".true")
    dfp.iter <- cbind(prop.true.iter, prop.pred.iter) %>% as.data.frame()
    Sys.sleep(system.sleep.sec)
    dfp.iter
  }))

  rownames(dfp) <- unique.sample.id.vector
  return(dfp)
}



dfp_tall_append_errors <- function(dfp.tall){
  # appends errors to a tall table
  # see also get_ypb_experiment_series()
  cn.tall <- colnames(dfp.tall)
  cn.tall <- cn.tall[!cn.tall %in% c("type", "sample.id", "abs.error.neuron")]
  ct.tall <- unique(gsub("\\..*", "", cn.tall)) 
  dfp.tall.new <- do.call(cbind, lapply(ct.tall, function(ct.iter){
    dfp.tall.iter <- dfp.tall[,grepl(paste0(ct.iter,"\\..*"), colnames(dfp.tall))]
    abs(dfp.tall.iter[,1]-dfp.tall.iter[,2])
  }))
  colnames(dfp.tall.new) <- paste0(ct.tall, ".abs.error")
  rownames(dfp.tall.new) <- rownames(dfp.tall)
  dfp.tall.new <- cbind(dfp.tall, dfp.tall.new)
  return(dfp.tall.new)
}

dfp_tall_by_celltype <- function(dfp.wide){
  # input: dfp.wide from get_ypb_experiment_series()
  # also append absolute errors 
  cn.wide <- colnames(dfp.wide)
  ct.wide <- unique(gsub("\\..*", "", cn.wide))
  dfp.tall.by.celltype <- do.call(rbind, lapply(ct.wide, function(ct.iter){
    dfp.wide.iter <- dfp.wide[,grepl(paste0(ct.iter,"\\..*"), cn.wide)]
    colnames(dfp.wide.iter) <- gsub(paste0(ct.iter, "\\."), "", colnames(dfp.wide.iter))
    dfp.wide.iter$celltype <- ct.iter
    dfp.wide.iter$sample.id <- rownames(dfp.wide.iter)
    dfp.wide.iter
  }))
  dfp.tall.by.celltype$abs.error.noscale <- abs(dfp.tall.by.celltype$true.noscale-dfp.tall.by.celltype$pred.noscale)
  dfp.tall.by.celltype$abs.error.withscale <- abs(dfp.tall.by.celltype$true.withscale-dfp.tall.by.celltype$pred.withscale)
  return(dfp.tall.by.celltype)
}
