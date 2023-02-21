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
num.iter <- 10
# set fract cells per iter
fract.cells.iter <- 20
# number of samples to select per iteration
num.sample.iter <- 3
# set the cell size factors for pseudobulk
S <- c("glial" = 3, "neuron" = 10)
# filepaths
rnf.data.fpath <- file.path("r-nf_deconvolution_donor-bias", "data")

# save paths
sce.fname <- "sce_all-samples.rda"
sce.fpath <- file.path(save.dpath, rnf.data.fpath, sce.fname)
# lindex
lindex.fname <- "lindex_inter.rda"
lindex.fpath <- file.path(save.dpath, rnf.data.fpath, lindex.fname)
# ypb
ypb.fname <- "ypb_inter.rda"
ypb.fpath <- file.path(save.dpath, rnf.data.fpath, ypb.fname)

# workflow table
wt.fname <- "workflow-table_inter-sample.csv"
wt.fpath <- file.path(save.dpath, rnf.data.fpath, wt.fname)

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
  scef <- sce[,vindex]
  
  # get true proportions vector
  tp <- as.data.frame(table(scef[[celltype.variable]]))
  tp.prop <- tp[,2]; names(tp.prop) <- tp[,1]
  tp.prop <- tp.prop/sum(tp[,2])
  
  # return results
  list(vindex = vindex, tp = tp.prop, ypb = ypb, samples = random.samples)
})
# save indices
save(lindex, file = lindex.fpath)

# get pseudobulk data
tp <- as.data.frame(table(sce[[celltype.variable]]))
tp.prop <- tp[,2]; names(tp.prop) <- tp[,1]
P <- tp.prop/sum(tp.prop)
Z <- do.call(cbind, lapply(unique.types, function(typei){
  rowMeans(assays(scef)[[assay.name]])
}))
ZS <- sweep(Z, 2, S, "*")
ypb <- t(t(P) %*% t(ZS))
# save pseudobulk
save(ypb, file = ypb.fpath)

#---------------------
# write workflow table
#---------------------
methodv <- c("nnls", "music")
wt <- do.call(rbind, lapply(methodv, function(methodi){
  wti <- data.frame(iterations_index = seq(num.iter))
  wti$method <- methodi
  wti$sample_id <- unlist(lapply(lindex, function(li){
    paste0(li$samples, collapse = ";")}))
  wti$celltype_variable <- celltype.variable
  wti$assay_name <- assay.name
  # manage filepaths
  cnamev <- c("sce_filepath", "bulk_filepath", "list_index_fpath")
  fpathv <- c(file.path("data", sce.fname),
              file.path("data", ypb.fname),
              file.path("data", lindex.fname))
  for(ii in seq(length(cnamev))){
    wti[,cnamev[ii]] <- paste0('"$launchDir/', fpathv[ii], '"')}
  wti
}))

# save
write.csv(wt, file = wt.fpath, row.names = F)









