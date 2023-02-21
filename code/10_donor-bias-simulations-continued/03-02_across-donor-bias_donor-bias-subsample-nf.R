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

# set number of iterations
num.iter.intra <- 10
num.iter.inter <- 10

# set fract cells per iter
fract.cells.iter <- 20

# filepaths
rnf.data.fpath <- file.path("r-nf_deconvolution_donor-bias", "data")

# intra data filepaths
# sce intra
sce.intra.fname <- "sce_largest-donor.rda"
sce.intra.fpath <- file.path(save.dpath, rnf.data.fpath, sce.intra.fname)
# lindex intra
lindex.fname <- "lindex_inter.rda"
lindex.fpath <- file.path(save.dpath, rnf.data.fpath, lindex.fname)

# workflow table
wt.fname <- "workflow-table_intra.csv"
wt.intra.fpath <- file.path(save.dpath, rnf.data.fpath, wt.fname)

#----------
# load data
#----------
fname <- paste0("list-scef_markers-k2-k3-k4_",proj.handle,".rda")
fpath <- file.path(load.dpath, fname)
lscef <- get(load(fpath))
sce <- lscef[[celltype.variable]]

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
num.donor.iter <- 3
cd <- colData(sce)
unique.samples <- unique(cd[,sample.variable])
unique.types <- unique(cd[,celltype.variable])
unique.types <- unique.types[order(unique.types)]
S <- c("glial" = 3, "neuron" = 10)
S <- S[order(match(names(S), unique.types))]
num.cellsv <- c("glial" = num.cells.glial, "neuron" = num.cells.neuron)
lindex <- lapply(seq(num.iter.inter), function(ii){
  set.seed(ii)
  random.samples <-  sample(unique.samples, num.donor.iter)
  filt <- cd[,sample.variable] %in% random.samples
  cdf <- cd[filt,]
  vindex <- which(colnames(sce) %in% unlist(lapply(unique.types, function(typei){
    sample(rownames(cdf[cdf[,celltype.variable]==typei,]), num.cellsv[typei])
  })))
  scef <- sce[,vindex]
  
  # get true proportions vector
  tp <- as.data.frame(table(scef[[celltype.variable]]))
  tp.prop <- tp[,2]; names(tp.prop) <- tp[,1]
  tp.prop <- tp.prop/sum(tp[,2])
  
  # get pseudobulk matrix
  P <- tp.prop
  # get signature matrix
  Z <- do.call(cbind, lapply(unique.types, function(typei){
    rowMeans(assays(scef)[[assay.name]])
  }))
  ZS <- sweep(Z, 2, S, "*")
  ypb <- t(t(P) %*% t(ZS))
  
  # return results
  list(vindex = vindex, tp = tp.prop, ypb = ypb, samples = random.samples)
})
# save indices
save(lindex, file = lindex.fpath)

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
wt <- data.frame(iterations_index = seq(num.iter.intra))
wt$method <- "nnls"
wt$sample_id <- largest.donor.id
wt$celltype_variable <- celltype.variable
wt$assay_name <- assay.name
# manage filepaths
cnamev <- c("sce_filepath", "true_proportions_filepath", 
            "bulk_filepath", "index_matrix_filepath")
fpathv <- c(file.path("data", sce.intra.fname),
            file.path("data", tp.intra.fname),
            file.path("data", ypb.intra.fname),
            file.path("data", mi.intra.fname))
for(ii in seq(length(cnamev))){wt[,cnamev[ii]] <- paste0('"$launchDir/', fpathv[ii], '"')}
# save
write.csv(wt, file = wt.intra.fpath, row.names = F)









