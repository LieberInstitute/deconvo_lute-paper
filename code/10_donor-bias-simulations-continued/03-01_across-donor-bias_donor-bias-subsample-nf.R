#!/usr/bin/env R

# Author: Sean Maden
#
# Run donor bias subsampling experiments.
#
# Note: run after script before 03-02 to use the same across-sample bulk 
# reference across the two experiments.
#

libv <- c("lute", "SummarizedExperiment", "SingleCellExperiment", 
          "ggplot2", "gridExtra")
sapply(libv, library, character.only = TRUE)

#------------------
# experiment params
#------------------
# get load path
code.dname <- "09_manuscript"
proj.dname <- "deconvo_method-paper"
load.dpath <- file.path(proj.dname, "outputs", code.dname)
# get save path
code.dname <- "10_donor-bias-simulations-continued"
proj.dname <- "deconvo_method-paper"
save.dpath <- file.path(proj.dname, "outputs", code.dname)
# get params for experiment
seed.num <- 0
celltype.variable <- "k2"
proj.handle <- "ro1-dlpfc"
save.fnstem <- paste0("inter-sample_", proj.handle)
group.variable <- "Sample"
assay.name <- "counts_adj"
methodv <- c("nnls", "music")
iterations <- 1000
fraction.cells <- 25
num.sample.iter <- 3
scale.factor <- c("glial" = 3, "neuron" = 10)
rnf.dname <- "r-nf_deconvolution_donor-bias"
base.path <- "data"
base.path <- file.path(save.dpath, rnf.dname, base.path)

#----------------------------------
# set up a new subsample experiment
#----------------------------------
# load data
fname <- paste0("list-scef_markers-k2-k3-k4_",proj.handle,".rda")
fpath <- file.path(load.dpath, fname)
lscef <- get(load(fpath))
sce <- lscef[[celltype.variable]]
rm(lscef)

# save new experiment data
lindex <- prepare_subsample_experiment(sce,
                                       scale.factor = scale.factor,
                                       iterations = iterations,
                                       groups.per.iteration = 3,
                                       method.vector = methodv,
                                       celltype.variable = celltype.variable,
                                       group.variable = group.variable,
                                       assay.name = assay.name,
                                       fraction.cells = fraction.cells,
                                       seed.num = seed.num,
                                       which.save = c("sce", "tp", "ypb", "li"),
                                       save.fnstem = save.fnstem,
                                       base.path = base.path,
                                       verbose = TRUE)

#---------------------
# manage workflow runs
#---------------------
# change wd
setwd(file.path(save.dpath, rnf.dname))
wt <- lindex$wt
# save full table
proj.handle <- "ro1-dlpfc"
save.fnstem <- paste0("inter-sample_", proj.handle)
wt.fnamei <- paste0("workflow-table-all_",save.fnstem,".csv")
write.csv(wt, file = file.path("data", wt.fnamei))

# from conda env with nextflow setup, run:
num.batch <- 500
indexv <- seq(1, nrow(wt), num.batch)
wt.fnamei <- paste0("workflow-table_",save.fnstem,".csv")
command.str <- paste0("bash ./sh/r-nf.sh -w ",
                      "workflow-table_inter-sample_ro1-dlpfc.csv")
lresults <- lapply(indexv, function(indexi){
  wt.path <- file.path("data", 
                       paste0("workflow-table-all_",save.fnstem,".csv"))
  wt <- read.csv(wt.path)
  wti <- wt[indexi:(indexi+num.batch-1),]
  write.csv(wti, file = file.path("data", wt.fnamei))
  system2(command.str)
  res.fname <- list.files()
  res.fname <- res.fname[grepl("results-table.*", res.fname)[1]]
  read.csv(res.fname)
})

for(indexi in indexv){
  
}


#---------------------
# running the workflow
#---------------------
## navigate to r-nf dir
# path = /deconvo_method-paper/outputs/10_donor-bias-simulations-continued/r-nf_deconvolution_donor-bias
# cd $path
#
## run the workflow:
# bash ./sh/r-nf.sh -w workflow-table_inter-sample_ro1-dlpfc.csv

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
