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

# list of z tables
lz.fname <- "lz-mr_expt-stransform_dlpfc-ro1.rda"
# sef.fname <- "sef-markers_stransform-expt_dlpfc-ro1.rda"
sef.fname <- "sef-markers_ct-treg_stransform-expt_dlpfc-ro1.rda"
script.fnamev <- c("z_methods.R", "z_transform.R", "make_example_data.R", 
                   "z_figures.R", "y_methods.R")

#-------
# params
#-------
celltype.treg.varname <- "celltype.treg"

#-----
# load
#-----
sef <- get(load(file.path(save.dpath, sef.fname)))

script.fnamev <- c("z_methods.R", "z_transform.R", "make_example_data.R", 
                   "z_figures.R", "y_methods.R")
# source key functions
for(scripti in script.fnamev){source(file.path(source.dpath, scripti))}

#----------------------
# get pseudobulk series
#----------------------
# set the cell weights
datv <- c(1,1,1,1)
# get the pseudobulked data
lpb <- get_lpb(scef = sef, datv = datv, get.results = TRUE,
               ctvarname = celltype.treg.varname, scale.range = 500:2000)
