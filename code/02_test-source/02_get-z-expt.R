#!/usr/bin/env R

#
# save the z expt
#

library(SingleCellExperiment)

#------
# paths
#------
proj.dpath <- "deconvo_method-paper"
sce.fpath <- "DLPFC_snRNAseq/processed-data/sce/sce_DLPFC.Rdata"
source.dpath <- file.path(proj.dpath, "source")
save.dpath <- file.path(proj.dpath, "outputs/02_test-source")
if(!dir.exists(save.dpath)){dir.create(save.dpath)}
lz.fname <- "lz_mr_dlpfc-ro1.rda"
script.fname <- "z_methods.R"

#----------
# load data
#----------
# single cell experiment
sce <- get(load(sce.fpath))
# source key functions
source(file.path(source.dpath, script.fname))

#-------------------
# get new zexpt data
#-------------------
# get detailed zexpt data
lz <- get_z_experiment(sce, method.markers = "mean_ratio", 
                       mr.assay = "logcounts", ngenes.byk = 25, 
                       type.varname = "cellType_broad_hc", 
                       summary.varname = "BrNum", k.summary.method = "mean",
                       z.summary.method = "mean", return.all = TRUE)

# save new expt data
save(lz, file = file.path(save.dpath, lz.fname))