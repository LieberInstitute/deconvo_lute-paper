#!/usr/bin/env R

library(here)
library(ggplot2)
library(ComplexHeatmap)

knitr::opts_chunk$set(echo = TRUE)
proj.dpath <- "deconvo_method-paper"
# source methods
source.dpath <- file.path("deconvo_method-paper", "source")
source.fnv <- c("pb_methods", "z_transformations")
for(fni in source.fnv){
  source(file.path(here(), source.dpath, paste0(fni, ".R")))}

#-----
# load
#-----
# get lz
lz.fname <- "lz_s-rescale-k4_dlpfc-ro1.rda"
lz.fpath <- file.path(proj.dpath, "outputs/03_test-stransform", lz.fname)
lz <- get(load(file.path(here(), lz.fpath)))
# get sef
sef.fname <- "sef-markers_ct-treg_stransform-expt_dlpfc-ro1.rda"
sef.fpath <- file.path(proj.dpath, "outputs/03_test-stransform", sef.fname)
sef <- get(load(file.path(here(), sef.fpath)))