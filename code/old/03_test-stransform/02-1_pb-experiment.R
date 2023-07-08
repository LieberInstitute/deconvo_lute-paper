#!/usr/bin/env R

#
# Running the stransformation experiment using more functions, updated pseudobulk data, etc.
# 

#------
# paths
#------
proj.dpath <- "deconvo_method-paper"
save.dpath <- file.path(proj.dpath, "outputs/03_test-stransform")
source.dpath <- file.path(proj.dpath, "source")
script.fnamev <- c("z_methods.R", "z_transform.R", 
                   "make_example_data.R", 
                   "z_figures.R", "y_methods.R")

#-----
# load
#-----
# source key functions
for(scripti in script.fnamev){
  source(file.path(source.dpath, scripti))}













