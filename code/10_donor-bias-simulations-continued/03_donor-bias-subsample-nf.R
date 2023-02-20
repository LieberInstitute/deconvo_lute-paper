#!/usr/bin/env R

# Author: Sean Maden
#
# Run donor bias subsampling experiments.
#
#

libv <- c("SummarizedExperiment", "SingleCellExperiment")
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
num.iter.intra <- 1000
num.iter.inter <- 1000

# set fract cells per iter
fract.cells.iter <- 20

# filepaths
rnf.data.fpath <- file.path("r-nf_deconvolution_donor-bias", "data")

# intra data filepaths
# sce intra
sce.intra.fname <- "sce_largest-donor.rda"
sce.intra.fpath <- file.path(save.dpath, rnf.data.fpath, sce.intra.fname)
# mindex intra
mi.fname <- "mindex_intra.rda"
mi.fpath <- file.path(save.dpath, rnf.data.fpath, mi.fname)
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
  datv <- sample(seq(num.neuron), num.neuron.select)
  datv <- c(datv, sample(seq(num.glial), num.glial.select))
}))
dim(mi) # [1] 1000 1398
save(mi, file = mi.fpath)

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

# write workflow table
wt <- data.frame(index = seq(1000))
wt$method <- "nnls"
wt$sample_id <- largest.donor.id
wt$celltype_variable <- celltype.variable
wt$assay_name <- assay.name
# manage filepaths
cnamev <- c("sce_filepath", "tp_filepath", "ypb_filepath")
fpathv <- c(file.path("data", sce.intra.fname),
            file.path("data", tp.intra.fname),
            file.path("data", ypb.intra.fname))
for(ci in seq(3)){wt[,cnamev[ci]] <- paste0('"$launchDir/', fpathv[ci], '"')}
# save
write.csv(wt, file = wt.intra.fpath, row.names = F)









