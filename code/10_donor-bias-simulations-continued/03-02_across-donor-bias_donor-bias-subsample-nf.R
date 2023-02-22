#!/usr/bin/env R

# Author: Sean Maden
#
# Run donor bias subsampling experiments.
#
#

libv <- c("SummarizedExperiment", "SingleCellExperiment", "ggplot2", "gridExtra")
sapply(libv, library, character.only = TRUE)

#------
# paths
#------
# get load path
code.dname <- "09_manuscript"
proj.dname <- "deconvo_method-paper"
load.dpath <- file.path(proj.dname, "outputs", code.dname)

# get save path
code.dname <- "10_donor-bias-simulations-continued"
proj.dname <- "deconvo_method-paper"
save.dpath <- file.path(proj.dname, "outputs", code.dname)

#-----------
# set params
#-----------
set.seed(0)
celltype.variable <- "k2"
proj.handle <- "ro1-dlpfc"
sample.variable <- "Sample"
assay.name <- "counts_adj"

# deconvolution methods to test
methodv <- c("nnls")

# set number of iterations
num.iter <- 400
# set fract cells per iter
fract.cells.iter <- 20
# number of samples to select per iteration
num.sample.iter <- 3
# set the cell size factors for pseudobulk
S <- c("glial" = 3, "neuron" = 10)
# filepaths
rnf.dname <- "r-nf_deconvolution_donor-bias"
rnf.data.fpath <- file.path(rnf.dname, "data")

# save paths
sce.fname <- "sce_all-samples.rda"
sce.fpath <- file.path(save.dpath, rnf.data.fpath, sce.fname)
# lindex
lindex.fname <- "lindex_inter.rda"
lindex.fpath <- file.path(save.dpath, rnf.data.fpath, lindex.fname)
# ypb
ypb.fname <- "ypb_inter.rda"
ypb.fpath <- file.path(save.dpath, rnf.data.fpath, ypb.fname)
# tp
tp.fname <- "true-proportions_inter.rda"
tp.fpath <- file.path(save.dpath, rnf.data.fpath, tp.fname)

# workflow table
wt.fname <- "workflow-table_inter-sample.csv"
wt.fpath <- file.path(save.dpath, rnf.data.fpath, wt.fname)

# results table
dfr.fname <- "results_table_inter.csv"
dfr.fpath <- file.path(save.dpath, rnf.dname, dfr.fname)

#----------
# load data
#----------
fname <- paste0("list-scef_markers-k2-k3-k4_",proj.handle,".rda")
fpath <- file.path(load.dpath, fname)
lscef <- get(load(fpath))
sce <- lscef[[celltype.variable]]

# save sce at r-nf.../data/
save(sce, file = sce.fpath)

#-------------------
# get random indices
#-------------------
# randomly take cells from 3 donors at a time
# get the exact numbers of cells of each type to sample
dft <- as.data.frame(table(sce[[celltype.variable]], 
                           sce[[sample.variable]]))
num.cells.glial <- round(min(dft[dft[,1]=="glial",3])*0.3, 0)
num.cells.neuron <- round(min(dft[dft[,1]=="neuron",3])*0.3, 0)
# get random cell indices as a list
cd <- colData(sce)
unique.samples <- unique(cd[,sample.variable])
unique.types <- unique(cd[,celltype.variable])
unique.types <- unique.types[order(unique.types)]
# order s cell size factors
S <- S[order(match(names(S), unique.types))]
num.cellsv <- c("glial" = num.cells.glial, "neuron" = num.cells.neuron)
lindex <- lapply(seq(num.iter), function(ii){
  set.seed(ii)
  # get filtered sce data as scef
  random.samples <-  sample(unique.samples, num.sample.iter)
  filt <- cd[,sample.variable] %in% random.samples
  cdf <- cd[filt,]
  vindex <- which(colnames(sce) %in% unlist(lapply(unique.types, function(typei){
    sample(rownames(cdf[cdf[,celltype.variable]==typei,]), num.cellsv[typei])
  })))
  # return results
  list(vindex = vindex, samples = random.samples)
})

# get pseudobulk data
tp <- as.data.frame(table(sce[[celltype.variable]]))
tp.prop <- tp[,2]; names(tp.prop) <- tp[,1]
P <- tp.prop/sum(tp.prop)
Z <- do.call(cbind, lapply(unique.types, function(typei){
  rowMeans(assays(sce)[[assay.name]])
}))
ZS <- sweep(Z, 2, S, "*")
ypb <- t(t(P) %*% t(ZS))

# save new data
# save indices
save(lindex, file = lindex.fpath)
# save pseudobulk
save(ypb, file = ypb.fpath)
# save tp
save(P, file = tp.fpath)

#---------------------
# write workflow table
#---------------------
wt <- do.call(rbind, lapply(methodv, function(methodi){
  wti <- data.frame(iterations_index = seq(num.iter))
  wti$method <- methodi
  wti$sample_id <- unlist(lapply(lindex, function(li){
    paste0(li$samples, collapse = ";")}))
  wti$celltype_variable <- celltype.variable
  wti$assay_name <- assay.name
  # manage filepaths
  cnamev <- c("sce_filepath", "bulk_filepath", "list_index_filepath",
              "true_proportions_filepath")
  fpathv <- c(file.path("data", sce.fname),
              file.path("data", ypb.fname),
              file.path("data", lindex.fname),
              file.path("data", tp.fname))
  for(ii in seq(length(cnamev))){
    wti[,cnamev[ii]] <- paste0('"$launchDir/', fpathv[ii], '"')}
  wti
}))

# save
write.csv(wt, file = wt.fpath, row.names = F)

#------------------------------
# load final results table data
#------------------------------
dfr <- read.csv(dfr.fpath)

#------------------------------------------------
# get final summary statistics from results table
#------------------------------------------------
methodv <- c(unique(dfr$deconvolution_method), "all")
funv <- c("median", "sd", "length")
metricv <- c("bias", "rmse.types")

# prepare results data
dfrs1 <- dfr
# add type label
cnv <- colnames(dfrs1)
unique.types <- unique(unlist(strsplit(dfr$type_labels, ";")))
unique.types <- unique.types[order(unique.types)]
for(typei in unique.types){
  cn.filt <- grepl(paste0("type", which(unique.types==typei)), cnv)
  colnames(dfrs1)[cn.filt] <- paste0(colnames(dfrs1)[cn.filt], ".", typei)
}
# get all method category
dfrs2 <- dfrs1; dfrs2$deconvolution_method <- "all"
dfrs3 <- rbind(dfrs1, dfrs2) # append all category
methods.vector <- dfrs3$deconvolution_method
lvar <- list(method = methods.vector)
unique.methods <- unique(methods.vector)
unique.methods <- unique.methods[order(unique.methods)]
# get new colnames for aggregate
cnv <- colnames(dfrs3)
grepl.str <- paste0(metricv, collapse = "|")
cnvf <- cnv[grepl(grepl.str, cnv)]

# get aggregate statistics
dfs <- do.call(cbind, lapply(funv, function(fi){
  dfai <- aggregate(dfrs3[,cnvf], lvar, FUN = fi); dfai <- dfai[,2:ncol(dfai)]
  fi.str <- fi; if(fi=="length"){fi.str <- "count"}
  colnames(dfai) <- paste0(fi.str, "_", colnames(dfai)); return(dfai)
}))
dfs$method <- unique.methods; cnv <- colnames(dfs)
dfs <- dfs[,c("method", cnv[!cnv=="method"])]

# save
fname <- "df-sstat-rnf_inter-donor-subsample_ro1-dlpfc.csv"
write.csv(dfs, file = file.path(save.dpath, fname), row.names = F)

#------------------------
# plot results table data
#------------------------
# scatter plots
# type bias by method
dfp <- data.frame(neuron = dfr$bias.type1, glial = dfr$bias.type2, 
                  method = dfr$deconvolution_method)
ggpt <- ggplot(dfp, aes(x = neuron, y = glial)) + theme_bw() +
  geom_point(alpha = 0.5, size = 2) + geom_abline(intercept = 0, slope = 1) + 
  xlim(-0.04, 0.04) + ylim(-0.04, 0.04)
ggpt + facet_wrap(~method)

# scatter plot -- bias
metric.plot <- title.str <- 'bias'
# type predictions by method
dfp1 <- dfr[dfr$deconvolution_method=="nnls",]
dfp2 <- dfr[dfr$deconvolution_method=="music",]
row1.orderv <- order(match(dfp1$iterations_index, dfp2$iterations_index))
dfp1 <- dfp1[row1.orderv,]
cond <- identical(dfp1$iterations_index, dfp2$iterations_index)
# get plot data
typev <- c(".type1", ".type2")
names(typev) <- unique(unlist(strsplit(dfp1$type_labels, ";")))
dfp <- do.call(rbind, lapply(c(".type1", ".type2"), function(typei){
  var.str <- paste0(metric.plot, typei)
  dfpi <- data.frame(nnls = dfp1[,var.str], music = dfp2[,var.str])
  dfpi$type <- names(typev[typev==typei]); dfpi
}))
# get plot object1
ggpt1 <- ggplot(dfp, aes(x = nnls, y = music)) + geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0) + theme_bw() + 
  geom_smooth() + ggtitle(title.str)
ggpt1 + facet_wrap(~type)
# get plot object2
ggpt2 <- ggplot(dfp, aes(x = nnls, y = music, color = type, shape = type)) + 
  geom_point(alpha = 0.5, size = 2) + geom_abline(slope = 1, intercept = 0) + 
  theme_bw() + geom_smooth() + ggtitle(title.str)
ggpt2

# scatter plot -- proportions
metric.plot <- title.str <- 'prop.pred'
# type predictions by method
dfp1 <- dfr[dfr$deconvolution_method=="nnls",]
dfp2 <- dfr[dfr$deconvolution_method=="music",]
row1.orderv <- order(match(dfp1$iterations_index, dfp2$iterations_index))
dfp1 <- dfp1[row1.orderv,]
cond <- identical(dfp1$iterations_index, dfp2$iterations_index)
# get plot data
typev <- c(".type1", ".type2")
names(typev) <- unique(unlist(strsplit(dfp1$type_labels, ";")))
dfp <- do.call(rbind, lapply(c(".type1", ".type2"), function(typei){
  var.str <- paste0(metric.plot, typei)
  dfpi <- data.frame(nnls = dfp1[,var.str], music = dfp2[,var.str])
  dfpi$type <- names(typev[typev==typei]); dfpi
}))
# get plot object1
ggpt1 <- ggplot(dfp, aes(x = nnls, y = music)) + geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0) + theme_bw() + 
  geom_smooth() + ggtitle(title.str)
ggpt1 + facet_wrap(~type)
# get plot object2
ggpt2 <- ggplot(dfp, aes(x = nnls, y = music, color = type, shape = type)) + 
  geom_point(alpha = 0.5, size = 2) + geom_abline(slope = 1, intercept = 0) + 
  theme_bw() + geom_smooth() + ggtitle(title.str)
ggpt2

# plots of rmse across types
title.str <- ""
variable.str <- ylab.str <- "rmse.types"
dfp <- dfr; dfp$value <- dfp[,variable.str]
# violin plots
ggvp1 <- ggplot(dfp, aes(x = deconvolution_method, y = value)) + 
  geom_violin(draw_quantiles = 0.5) + theme_bw() +
  ggtitle(title.str) + ylab(ylab.str)
# jitter plots
ggjt1 <- ggplot(dfp, aes(x = deconvolution_method, y = value)) + 
  geom_jitter(alpha = 0.5) + geom_boxplot(alpha = 0, color = "purple", lwd = 1) + 
  theme_bw() + facet_zoom(ylim = c(0, 1e-16)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle(title.str) + ylab(ylab.str)

# plot of rmse by type
get_rmse_type <- function(prop.true, prop.pred){
  unlist(lapply(seq(length(prop.true)), function(ii){
    sqrt(mean((prop.pred[ii]-prop.true[ii])^2))
  }))
}
typev <- unique(unlist(strsplit(dfr$type_labels, ";")))
dfp <- do.call(rbind, lapply(c("nnls", "music"), function(methodi){
  dfri <- dfr[dfr$deconvolution_method==methodi,]
  dfpi <- do.call(cbind, lapply(typev, function(typei){
    type.index <- which(typev==typei)
    varname <- paste0("type", type.index)
    varname.pred <- paste0("prop.pred.", varname)
    varname.true <- paste0("prop.true.", varname)
    get_rmse_type(dfri[,varname.true], dfri[,varname.pred])
  }))
  dfpi <- as.data.frame(dfpi)
  colnames(dfpi) <- typev
  dfpi$method = methodi; dfpi
}))

ylab.str <- "RMSE"
lgg <- lapply(typev, function(typei){
  variable.name <- title.str <- typei
  dfp$value <- dfp[,variable.name]
  # violin plots 
  ggvp1 <- ggplot(dfp, aes(x = method, y = value)) + 
    geom_violin(draw_quantiles = 0.5) + theme_bw() +
    ggtitle(title.str) + ylab(ylab.str)
  # jitter plots
  ggjt1 <- ggplot(dfp, aes(x = method, y = value)) + 
    geom_jitter(alpha = 0.5) + 
    geom_boxplot(alpha = 0, color = "purple", lwd = 1) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(title.str) + ylab(ylab.str)
  list(violin = ggvp1, jitter = ggjt1)
})
names(lgg) <- typev
lgg$neuron$violin
lgg$glial$violin
lgg$neuron$jitter
lgg$glial$jitter
