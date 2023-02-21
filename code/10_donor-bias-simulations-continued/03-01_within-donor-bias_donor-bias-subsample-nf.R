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

#--------------
# load data
#--------------
proj.handle <- "ro1-dlpfc"
fname <- paste0("list-scef_markers-k2-k3-k4_",proj.handle,".rda")
fpath <- file.path(load.dpath, fname)
lscef <- get(load(fpath))
sce <- lscef[["k2"]]

#---------------
# set params
#---------------
set.seed(0)
celltype.variable <- "k2"
donor.variable <- "Sample"
assay.name <- "counts_adj"

# set number of iterations
num.iter.intra <- num.iter.inter <- 100

# set fract cells per iter
fract.cells.iter <- 20

# filepaths
rnf.fpath <- "r-nf_deconvolution_donor-bias"
rnf.data.fpath <- file.path(rnf.fpath, "data")

# intra data filepaths
# sce intra
sce.intra.fname <- "sce_largest-donor.rda"
sce.intra.fpath <- file.path(save.dpath, rnf.data.fpath, sce.intra.fname)
# mindex intra
mi.intra.fname <- "mindex_intra.rda"
mi.intra.fpath <- file.path(save.dpath, rnf.data.fpath, mi.intra.fname)
# true proportions
tp.intra.fname <- "true-proportions_intra.rda"
tp.intra.fpath <- file.path(save.dpath, rnf.data.fpath, tp.intra.fname)
# pseudobulk
ypb.intra.fname <- "ypb_intra.rda"
ypb.intra.fpath <- file.path(save.dpath, rnf.data.fpath, ypb.intra.fname)
# workflow table
wt.fname <- "workflow-table_intra.csv"
wt.intra.fpath <- file.path(save.dpath, rnf.data.fpath, wt.fname)

#------------------------------------
# subsample, save within-donor biases
#------------------------------------
# get largest donor
dft <- as.data.frame(table(sce[[donor.variable]]))
largest.donor.id <- dft[dft[,2]==max(dft[,2]),1]
scef <- sce[,sce[[donor.variable]]==largest.donor.id]
# save
save(scef, file = sce.intra.fpath)

# save indices
celltype.variable.vector <- scef[[celltype.variable]]
num.neuron <- length(which(celltype.variable.vector == "neuron"))
num.glial <- length(which(celltype.variable.vector == "glial"))
num.neuron.select <- round(num.neuron * fract.cells.iter/100, 0)
num.glial.select <- round(num.neuron * fract.cells.iter/100, 0)
mi <- do.call(rbind, lapply(seq(num.iter.intra), function(ii){
  set.seed(ii)
  which.neuron.iter <- which(scef[[celltype.variable]] == "neuron")
  which.neuron.iter <- which.neuron.iter[sample(seq(num.neuron), num.neuron.select)]
  which.glial.iter <- which(scef[[celltype.variable]] == "glial")
  which.glial.iter <- which.glial.iter[sample(seq(num.glial), num.glial.select)]
  c(which.neuron.iter, which.glial.iter)
}))
colnames(mi) <- c(rep("neuron", num.neuron.select), rep("glial", num.glial.select))
dim(mi) # [1] 1000 1398
save(mi, file = mi.intra.fpath)

# true proportions
dft <- as.data.frame(table(scef[[celltype.variable]]))
tp <- dft[,2]/sum(dft[,2])
names(tp) <- dft[,1]
save(tp, file = tp.intra.fpath)

# pseudobulk
# true proportions 
P <- tp
# get signature matrix
unique.types <- unique(scef[[celltype.variable]])
Z <- do.call(cbind, lapply(unique.types, function(typei){
  rowMeans(assays(scef)[[assay.name]])
}))
S <- c(10, 3)
ZS <- sweep(Z, 2, S, "*")
ypb <- t(t(P) %*% t(ZS))
# save 
save(ypb, file = ypb.intra.fpath)

#---------------------
# write workflow table
#---------------------
# note: make table around number of methods to try
methodv <- c("nnls", "music")
wt <- do.call(rbind, lapply(methodv, function(methodi){
  wti <- data.frame(iterations_index = seq(num.iter.intra))
  wti$method <- methodi
  wti$sample_id <- largest.donor.id
  wti$celltype_variable <- celltype.variable
  wti$assay_name <- assay.name
  # manage filepaths
  cnamev <- c("sce_filepath", "true_proportions_filepath", 
              "bulk_filepath", "index_matrix_filepath")
  fpathv <- c(file.path("data", sce.intra.fname),
              file.path("data", tp.intra.fname),
              file.path("data", ypb.intra.fname),
              file.path("data", mi.intra.fname))
  for(ii in seq(length(cnamev))){
    wti[,cnamev[ii]] <- paste0('"$launchDir/', fpathv[ii], '"')}
  wti
}))

# save
write.csv(wt, file = wt.intra.fpath, row.names = F)

#----------------------
# analyze results table
#----------------------
results.fpath <- "results_table_intra.csv"
dfr <- read.csv(file.path(save.dpath, rnf.fpath, results.fpath))

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

# plots of bias, rmse distributions by method
# violin plots
ggvp1 <- ggplot(dfr, aes(x = deconvolution_method, y = rmse.types)) + 
  geom_violin(draw_quantiles = 0.5) + theme_bw()
# jitter plots
ggjt1 <- ggplot(dfr, aes(x = deconvolution_method, y = rmse.types)) + 
  geom_jitter(alpha = 0.5) + stat_summary() + theme_bw()

#-----------------------------
# get final summary statistics
#-----------------------------

