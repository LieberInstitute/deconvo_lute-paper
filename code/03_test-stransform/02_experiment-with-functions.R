#!/usr/bin/env R

#
# Running the stransformation experiment using more functions, updated pseudobulk data, etc.
# 

#------
# paths
#------
proj.dpath <- "deconvo_method-paper"
source.dpath <- file.path(proj.dpath, "source")
save.dpath <- file.path(proj.dpath, "outputs/03_test-stransform")

lz.fname <- "lz-mr_expt-stransform_dlpfc-ro1.rda"

# sce.fpath <- "DLPFC_snRNAseq/processed-data/sce/sce_DLPFC.Rdata"# 
sef.fname <- "sef-markers_stransform-expt_dlpfc-ro1.rda"

script.fnamev <- c("z_methods.R", "z_transform.R", 
                   "make_example_data.R", "z_figures.R", 
                   "y_methods.R")

# saved objects
# new z datasets
lz.fname <- "lz_s-rescale-k4_dlpfc-ro1.rda"
# new pseudobulked data
lpb.fname <- "lpseudobulk_stransform-expt_dlpfc-ro1.rda"

#-----
# load
#-----
sef <- get(load(file.path(save.dpath, sef.fname)))

