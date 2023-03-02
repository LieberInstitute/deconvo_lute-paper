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
iterations <- 100
# fraction.cells <- 25
num.sample.iter <- 3
scale.factor <- c("glial" = 3, "neuron" = 10)
rnf.dname <- "r-nf_deconvolution_donor-bias"
base.path <- "data"
base.path <- file.path(save.dpath, rnf.dname, base.path)
# save data
which.save = c("li", "sce")
save.names = list(sce.name = "sce.rda",
                  li.name = "lindex.rda")

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
# get main starting workflow table
wt <- lsub$wt
# save full table
proj.handle <- "ro1-dlpfc"
save.fnstem <- paste0("intra-sample_", proj.handle)
wt.fnamei <- paste0("workflow-table-all_",save.fnstem,".csv")
wt.fpath <- file.path(save.dpath, rnf.dname, "data", wt.fnamei)
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
  wti.fpath <- file.path(save.dpath, rnf.dname, "data", wti.fname.iter)
  write.csv(wti, file = wti.fpath, row.names = F)
  message("finished with iter ", iter.num)
}

#---------------------------------------
# get the point estimate of bias by type -- no cell size factors
#---------------------------------------
base.path <- "data"
base.path <- file.path(save.dpath, rnf.dname, base.path)

# get true proportions
P <- prop.table(table(scef[[celltype.variable]]))
# get signature matrix
unique.types <- unique(scef[["k2"]])
unique.types <- unique.types[order(unique.types)]
mexpr <- assays(scef)[["counts_adj"]]
Z <- do.call(cbind, lapply(unique.types, function(typei){
  filt <- scef[["k2"]]==typei; rowMeans(mexpr[,filt])}))
colnames(Z) <- unique.types
Z <- as.data.frame(Z)

# get ypb -- no cell scale
ypb <- t(t(P) %*% t(Z))
# get pred
arg <- list()
ld2 <- run_deconvolution(Z = as.matrix(Z), Y = as.matrix(ypb), method = "music", arguments = arg)
ld3 <- run_deconvolution(Z = as.matrix(Z), Y = as.matrix(ypb), method = "nnls", arguments = arg)
ld1$predictions
ld2$predictions
ld3$predictions

# get ypb -- with cell size scale on Y but not Z
ZS <- sweep(Z, 2, S, "*")
ypb <- t(t(P) %*% t(ZS))
# get pred
arg <- list()
ld1 <- run_deconvolution(Z = Z, Y = ypb, method = "epic", arguments = arg)
ld2 <- run_deconvolution(Z = as.matrix(Z), Y = as.matrix(ypb), method = "music", arguments = arg)
ld3 <- run_deconvolution(Z = as.matrix(Z), Y = as.matrix(ypb), method = "nnls", arguments = arg)
ld1$predictions
ld2$predictions
ld3$predictions

# get ypb -- with cell size scale on both Y and Z
ZS <- sweep(Z, 2, S, "*")
ypb <- t(t(P) %*% t(ZS))
# Z <- sweep(Z, 2, S, "*")
# get pred
arg <- list()
Zvar <- do.call(cbind, lapply(unique.types, function(typei){
  rowVars(mexpr[,scef[["k2"]]==typei])
}))
rownames(Zvar) <- rownames(Z); colnames(Zvar) <- colnames(Z)
reference <- list(refProfiles = as.matrix(ZS), sigGenes = rownames(ZS),
                  refProfiles.var = as.matrix(Zvar))
ld1 <- EPIC(ypb, reference = reference, mRNA_cell = S)
ld1$mRNAProportions

ld1 <- run_deconvolution(Z = Z, Y = ypb, method = "epic", arguments = arg)
ld2 <- run_deconvolution(Z = as.matrix(Z), Y = as.matrix(ypb), 
                         method = "music", arguments = arg)
ld3 <- run_deconvolution(Z = as.matrix(Z), Y = as.matrix(ypb), 
                         method = "nnls", arguments = arg)
ld1$predictions
ld2$predictions
ld3$predictions

#--------------------------
# test subsample parameters
#--------------------------
# define main function
get_bias <- function(ypb, sce, S = NULL, num.cells = 100, 
                     num.iter = 10, method = "music", seed.num = 0,
                     proportions = c("glial" = 0.3, "neuron" = 0.7)){
  set.seed(seed.num)
  df.prop <- do.call(rbind, lapply(seq(num.iter), function(ii){
    message("iteration: ", ii) # ; set.seed(ii)
    Znew <- get_random_z(proportions, num.cells, sce, seed.num = ii)
    if(!is(S, "NULL")){Znew <- sweep(Znew, 2, S, "*")}
    ld2 <- suppressMessages(run_deconvolution(Z = as.matrix(Znew), 
                                              Y = as.matrix(ypb), 
                                              method = method))
    ld2$predictions
  }))
  df.bias <- sweep(df.prop, 1, proportions, "-")
  return(df.bias)
}

get_random_z <- function(proportions, num.cells, sce, celltype.variable = "k2", 
                         assay.name = "counts_adj", seed.num = 0){
  set.seed(seed.num)
  unique.types = names(proportions)[order(names(proportions))]
  mexpr <- assays(sce)[[assay.name]]
  cell.num <- get_cell_quantities(proportions, num.cells)
  Znew <- do.call(cbind, lapply(unique.types, function(typei){
    sce.indexv <- which(sce[[celltype.variable]]==typei)
    type.indexv <- sample(seq(sce.indexv), cell.num[typei])
    mef <- mexpr[,sce.indexv[type.indexv]]
    rowMeans(mef)
  }))
  rownames(Znew) <- rownames(sce)
  colnames(Znew) <- unique.types
  return(Znew)
}

# get ypb
S <- c("glial" = 3, "neuron" = 10)
# get true proportions
P <- prop.table(table(scef[[celltype.variable]]))
# get signature matrix
unique.types <- unique(scef[["k2"]])
unique.types <- unique.types[order(unique.types)]
mexpr <- assays(scef)[["counts_adj"]]
Z <- do.call(cbind, lapply(unique.types, function(typei){
  filt <- scef[["k2"]]==typei; rowMeans(mexpr[,filt])}))
colnames(Z) <- unique.types
Z <- as.data.frame(Z)
ZS <- sweep(Z, 2, S, "*")
ypb <- t(t(P) %*% t(ZS))# get ypb -- no cell scale

# test bias function
df.bias.exe <- get_bias(ypb, scef, S = S)
# plot
df.bias.exe <- as.data.frame(df.bias.exe)
colnames(df.bias.exe) <- unique.types
ggplot(df.bias.exe, aes(x = glial, y = neuron)) + 
  geom_abline(slope = 1, intercept = 0) + geom_point(alpha = 0.4)

# get data
num.iter <- 2

P <- c("glial" = 0.3, "neuron" = 0.7)
num.cells.vector <- seq(20, 1000, 200)
dfs <- do.call(rbind, lapply(num.cells.vector, function(num.cells){
  dat <- suppressMessages(get_bias(ypb, S = S, num.cells = num.cells, 
                            num.iter = num.iter, method = "nnls",
                            proportions = P))
}))
dfs <- as.data.frame(dfs)
dfs$s.factor <- TRUE
dfs2 <- do.call(rbind, lapply(num.cells.vector, function(num.cells){
  suppressMessages(get_bias(ypb, S = NULL, num.cells = num.cells, 
                            num.iter = num.iter, method = "nnls",
                            proportions = P))
}))
dfs2 <- as.data.frame(dfs2)
dfs2$s.factor <- FALSE
# get plot data
dfp <- rbind(dfs, dfs2)
colnames(dfp) <- c("bias.glial", "bias.neuron", "s.factor")
dfp$num.cells <- rep(rep(num.cells.vector, each = num.iter), 2)
dfp$abs.bias.glial <- abs(dfp$bias.glial)
# save dfp
#dfp.fname <- "dfp-subsample-param_intra-sample-bias.rda"
#save(dfp, file = file.path(save.dpath, dfp.fname))

ggplot(dfp, aes(x = num.cells, y = bias.glial, group = s.factor, color = s.factor)) +
  geom_hline(yintercept = 0, color = "black") +
  geom_point() + geom_smooth()

ggplot(dfp, aes(x = num.cells, y = abs.bias.glial, group = s.factor, color = s.factor)) +
  geom_hline(yintercept = 0, color = "black") +
  geom_point() + geom_smooth()

ggplot(dfp, aes(x = num.cells, y = bias.glial, group = s.factor, color = s.factor)) +
  geom_hline(yintercept = 0, color = "black") +
  geom_point() + geom_smooth() + facet_grid(~s.factor)

ggplot(dfp, aes(x = s.factor, y = bias.glial)) +
  geom_hline(yintercept = 0, color = "black") +
  geom_violin()

ggplot(dfp, aes(x = s.factor, y = bias.glial)) + geom_jitter() + 
  geom_hline(yintercept = 0, color = "black") +
  geom_boxplot(alpha = 0)

dfp$num.cells <- as.factor(dfp$num.cells)
ggplot(dfp, aes(x = num.cells, y = bias.glial, color = s.factor)) + 
  geom_hline(yintercept = 0, color = "black") +
  geom_jitter() + geom_boxplot(alpha = 0)

ggplot(dfp, aes(x = bias.glial, y = bias.neuron, color = s.factor)) + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  geom_point()

# try with random cell subset
set.seed(0)
num.cells <- table(scef[["k2"]])
num.glial <- num.cells["glial"]
num.neuron <- num.cells["neuron"]
dbias <- do.call(rbind, lapply(seq(10), function(ii){
  message(ii); set.seed(ii)
  # get random cells
  cell.num <- get_cell_quantities(c("glial" = 0.3, "neuron" = 0.7), 400)
  glial.indexv <- sample(seq(num.glial), cell.num["glial"])
  glial.indexv <- which(scef[["k2"]]=="glial")[glial.indexv]
  glial.expr <- rowMeans(mexpr[,glial.indexv])
  neuron.indexv <- sample(seq(num.neuron), cell.num["neuron"])
  neuron.indexv <- which(scef[["k2"]]=="neuron")[neuron.indexv]
  neuron.expr <- rowMeans(mexpr[,neuron.indexv])
  # get decon results
  Znew <- cbind(glial.expr, neuron.expr)
  ZSnew <- sweep(Znew, 2, S, "*")
  print(head(ZSnew))
  ld2 <- suppressMessages(run_deconvolution(Z = as.matrix(ZSnew), 
                                            Y = as.matrix(ypb), 
                           method = "nnls", arguments = arg))
  ld2$predictions
}))

dbias2 <- sweep(dbias, 1, P, "-")
dbias2
median(abs(dbias2[,1]))
median(abs(dbias2[,2]))

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
jpeg("ggvp-comp_intra-sample-bias.jpg", width = 5, height = 4,
     units = "in", res = 400)
grid.arrange(lgg[[1]], lgg[[2]], lgg[[3]], lgg[[4]], lgg[[5]],
             layout_matrix = matrix(c(1,2,3,4,5), nrow = 2))
dev.off()
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
