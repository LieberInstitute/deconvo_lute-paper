#!/usr/bin/env R

#
# Make the filtered summarized experiment, non-DelayedArray.
# 

library(SingleCellExperiment)
library(SummarizedExperiment)

#------
# paths
#------
proj.dpath <- "deconvo_method-paper"
save.dpath <- file.path(proj.dpath, "outputs/03_test-stransform")
sce.fpath <- "DLPFC_snRNAseq/processed-data/sce/sce_DLPFC.Rdata"
lz.fname <- "lz_s-rescale-k4_dlpfc-ro1.rda"
sef.fname <- "sef-markers_ct-treg_stransform-expt_dlpfc-ro1.rda"

#-----
# load
#-----
# single cell experiment
sce <- get(load(sce.fpath))
# gene markers
lz <- get(load(file.path(save.dpath, lz.fname)))
genemarkerv <- rownames(lz[[2]])

#------
# param
#------
celltype.varname <- "cellType_broad_hc"
celltype.treg.varname <- "celltype.treg"

#--------------
# remake as sef
#--------------
scef <- sce[genemarkerv,]
lct <- list(counts = as.matrix(counts(scef)))
sef <- SummarizedExperiment(assays = lct)
colData(sef) <- colData(sce)

#---------------------
# make tregs cell type
#---------------------
varv <- as.character(sce[[celltype.varname]])
varv[!varv %in% c("Excit", "Inhib", "Oligo")] <- "other"
sef[[celltype.treg.varname]] <- varv

#-----
# save
#-----
save(sef, file = file.path(save.dpath, sef.fname))
