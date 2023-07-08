#!/usr/bin/env R

# Author: Sean Maden 
#
# Prepping multiregion brain dataset for analysis.

libv <- c("SingleCellExperiment", "SummarizedExperiment")
sapply(libv, library, character.only = TRUE)

# manage paths
code.dname <- "08_lute-simulations"
proj.dname <- "deconvo_method-paper"
save.dpath <- file.path(proj.dname, "outputs", code.dname)
source.dpath <- file.path(proj.dname, "source")
script.fname <- "helperfun_tran-et-al.R"
source(file.path(source.dpath, script.fname))

#--------------
# prep sce data
#--------------

# get data
sce.fname <- "sce-mrb_dlpfc.rda"
sce.fpath <- file.path(save.dpath, sce.fname)
sce <- get_sce_mrb()
sce <- prep_sce()
save(sce, file = sce.fpath)

# check attributes
names(assays(sce))
# [1] "counts"    "logcounts"
table(sce[["k2"]])
# neuron  other 
# 3968   7234

#-------------------
# get marker objects
#-------------------
# get markers of cell types
mr <- get_mean_ratio2(sce, "k2", "logcounts")
mr.fname <- "markers-k2_db-mr2_sce-dlpfc-mrb.rda"
mr.fpath <- file.path(save.dpath, mr.fname)
save(mr, file = mr.fpath)

# get top markers -- 20 genes/type
mr.fname <- "markers-k2_db-mr2-t20_sce-dlpfc-mrb.rda"
mr.fpath <- file.path(save.dpath, mr.fname)
type.variable <- "k2"
ngenes.per.k <- 20
typev <- unique(sce[["k2"]])
ma.top <- do.call(rbind, lapply(typev, function(typei){
  mr %>% filter(cellType.target == typei) %>% arrange(rank_ratio) %>%
    top_n(n = ngenes.per.k)
}))
# get scef
scef <- sce[ma.top$gene,]
scef.fname <- "scef-mrb_gmr2-k2-t20_dlpfc.rda"
save(scef, file = file.path(save.dpath, scef.fname))

# get markers of donors
mrd <- get_mean_ratio2(sce, "donor", "logcounts")
mrd.fname <- "markers-donor_db-mr2_sce-dlpfc-mrb.rda"
mrd.fpath <- file.path(save.dpath, mrd.fname)
save(mrd, file = mrd.fpath)




