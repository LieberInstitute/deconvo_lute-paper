#!/usr/bin/env R

# Author: Sean Maden
#
# Run donor bias subsampling experiments.
# 
# Note: run after script 03-01 to use the same across-sample bulk reference.
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
save.fnstem <- paste0("intra-sample_", proj.handle)
group.variable <- "Sample"
assay.name <- "counts_adj"
methodv <- c("nnls", "music", "epic", "deconrnaseq")
iterations <- 2000
# fraction.cells <- 25
num.sample.iter <- 3
scale.factor <- c("glial" = 3, "neuron" = 10)
rnf.dname <- "r-nf_deconvolution_intra-sample-bias"
base.path <- "data"
base.path <- file.path(save.dpath, rnf.dname, base.path)
# save data
which.save = c("li", "sce")
save.names = list(sce.name = "sce.rda", li.name = "lindex.rda")

#----------------------------------
# set up a new subsample experiment
#----------------------------------
# load data
fname <- paste0("list-scef_markers-k2-k3-k4_",proj.handle,".rda")
fpath <- file.path(load.dpath, fname)
lscef <- get(load(fpath))
sce <- lscef[[celltype.variable]]
rm(lscef)

# get largest donor
dft <- as.data.frame(table(sce[[group.variable]]))
largest.donor.id <- dft[dft[,2]==max(dft[,2]),1]
scef <- sce[,sce[[group.variable]]==largest.donor.id]

# set cell count minimum
count.min <- 200
# get proportions of types in scef
proportions <- prop.table(table(scef[["k2"]]))
proportions
# glial    neuron 
# 0.2489252 0.7510748

# save new experiment data
lsub <- prepare_subsample_experiment(scef,
                                     scale.factor = scale.factor,
                                     iterations = iterations,
                                     groups.per.iteration = 1,
                                     method.vector = methodv,
                                     celltype.variable = celltype.variable,
                                     group.variable = group.variable,
                                     assay.name = assay.name,
                                     cell.proportions = proportions,
                                     count.minimum = count.min,
                                     seed.num = seed.num,
                                     which.save = c("sce", "tp", "ypb", "li"),
                                     save.fnstem = save.fnstem,
                                     base.path = "data",
                                     verbose = TRUE)

#---------------------
# manage workflow runs
#---------------------
base.path <- "data" # file.path(save.dpath, rnf.dname, "data")
# get main starting workflow table
wt <- lsub$wt
# save full table
proj.handle <- "ro1-dlpfc"
save.fnstem <- paste0("intra-sample_", proj.handle)
wt.fnamei <- paste0("workflow-table-all_",save.fnstem,".csv")
wt.fpath <- file.path(base.path, wt.fnamei)
write.csv(wt, file = wt.fpath)

# save table iterations
num.batch <- 200
indexv <- seq(1, nrow(wt), num.batch)
wt.fnamei <- paste0("workflow-table_",save.fnstem,".csv")
for(iter.num in seq(length(indexv))){
  index.start <- indexv[iter.num]
  index.end <- (index.start+num.batch-1)
  if(index.end > nrow(wt)){index.end <- nrow(wt)}
  wti <- wt[index.start:index.end,]
  wti.fname.iter <- wt.fnamei <- paste0("workflow-table-iter", 
                                        iter.num, "_", save.fnstem,".csv")
  wti.fpath <- file.path(base.path, wti.fname.iter)
  write.csv(wti, file = wti.fpath, row.names = F)
  message("finished with iter ", iter.num)
}

#-----------------
# run the workflow
#-----------------
# navigate to:
# /deconvo_method-paper/outputs/10_donor-bias-simulations-continued/r-nf_deconvolution_donor-bias
#
# run:
# bash ./sh/r-nf_run_wt-iter.sh

#-----------------------
# analyze results table
#-----------------------
results.filt <- "results_table_.*"
data.dpath <- file.path(save.dpath, rnf.dname, "data")
lfv <- list.files(data.dpath)
rt.fname <- lfv[grepl(results.filt, lfv)]
rt.fpath <- file.path(save.dpath, rnf.dname, "data", rt.fname)
rt <- read.csv(rt.fpath)


library(ggplot2)
library(gridExtra)
library(GGally)

# plot pairs
method.varname <- "deconvolution_method"
# varname <- "prop.pred.type1"
# filter duplicate entries
rtf <- rt
rtf$label <- paste0(rtf$iterations_index, ";", rtf$deconvolution_method)
rtf <- rtf[!duplicated(rtf$label),]
methodv <- unique(rtf[,method.varname])
# iterate on plot variables
varv <- c("prop.pred.type1", "prop.pred.type2", "bias.type1", "bias.type2", "rmse.types")
lgg <- lapply(varv, function(varname){
  dfp <- do.call(cbind, lapply(methodv, function(methodi){
    filt <- rtf[,method.varname]==methodi; rti <- rtf[filt,]
    rti <- rti[order(rti$iterations_index),]; rti[,varname]
  }))
  dfp <- as.data.frame(dfp); colnames(dfp) <- methodv
  ggpairs(dfp) + ggtitle(varname)
})
names(lgg) <- varv

lgg$prop.pred.type1
lgg$prop.pred.type2
lgg$bias.type1
lgg$bias.type2
lgg$rmse.types

# violin plots
#
lgg <- lapply(varv, function(varname){
  dfp <- rtf; dfp$value <- dfp[,varname]
  ggplot(dfp, aes(x = deconvolution_method, y = value)) +
    geom_violin(draw_quantiles = 0.5) + ggtitle(varname)
})
# bias plots
plot1 <- lgg[[3]] + theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  geom_hline(yintercept = 0, alpha = 0.2)
plot2 <- lgg[[4]] + theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0, alpha = 0.2)
jpeg("ggvp-comp-bias_intra-sample-bias.jpg", 
     width = 6, height = 5, units = "in", res = 400)
grid.arrange(plot1, plot2, nrow = 2,
             layout_matrix = matrix(c(1,1,1,2,2,2,2), ncol = 1),
             bottom = "Deconvolution method")
dev.off()

# jitter plots
lgg <- lapply(varv, function(varname){
  dfp <- rtf; dfp$value <- dfp[,varname]
  ggplot(dfp, aes(x = deconvolution_method, y = value)) +
    geom_jitter(alpha = 0.5) + ggtitle(varname) +
    geom_boxplot(alpha = 0, color = "cyan", size = 1)
})
# bias plots
plot1 <- lgg[[3]] + theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  geom_hline(yintercept = 0, alpha = 1, color = "red")
plot2 <- lgg[[4]] + theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0, alpha = 1, color = "red")
jpeg("ggjitterbox-comp-bias_intra-sample-bias.jpg", 
     width = 6, height = 5, units = "in", res = 400)
grid.arrange(plot1, plot2, nrow = 2,
             layout_matrix = matrix(c(1,1,1,2,2,2,2), ncol = 1),
             bottom = "Deconvolution method")
dev.off()

# scatter plot bias
# no facet
ggpt <- ggplot(rtf, aes(x = bias.type1, y = bias.type2, 
                        group = deconvolution_method, 
                        color = deconvolution_method)) +
  geom_point(alpha = 0.2) + geom_abline(slope = 1, intercept = 0) +
  theme_bw() + 
  geom_hline(yintercept = 0, alpha = 0.5, color = "gray") +
  geom_vline(xintercept = 0, alpha = 0.5, color = "gray")
# facet
ggpt <- ggplot(rtf, aes(x = bias.type1, y = bias.type2)) +
  geom_point(alpha = 0.2) + geom_abline(slope = 1, intercept = 0) +
  theme_bw() + 
  geom_hline(yintercept = 0, alpha = 0.5, color = "gray") +
  geom_vline(xintercept = 0, alpha = 0.5, color = "gray")
ggpt <- ggpt + facet_wrap(~deconvolution_method)
jpeg('ggpt-facet-method_intra-sample-bias.jpg', 
     width = 6, height = 6, 
     units = "in", res = 400)
ggpt
dev.off()
