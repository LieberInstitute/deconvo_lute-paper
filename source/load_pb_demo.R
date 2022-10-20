#!/usr/bin/env R

library(here)
# library(ggplot2)
knitr::opts_chunk$set(echo = TRUE)
message("here: ", here())
proj.dpath <- "deconvo_method-paper"
# source methods
source.dpath <- file.path("deconvo_method-paper",
                          "source")
source.fnv <- c("pb_methods", 
                "z_transformations")
for(fni in source.fnv){
  source(
    file.path(here(), 
              source.dpath, 
              paste0(fni, ".R")))}
lz.fname <- "lz_s-rescale-k4_dlpfc-ro1.rda"
lz.fpath <- file.path(proj.dpath, 
                      "outputs/03_test-stransform", 
                      lz.fname)
lz <- get(load(file.path(here(), lz.fpath)))
#dft <- data.frame(cell_type = c("Excit", "Inhib", "Oligo", "other"),
#                        mean = c(6, 8, 2, 2), 
#                        sd = c(2, 3, 1, 1))
#meanv <- c(6, 2, 2, 8)
#sdv <- c(2, 1, 1, 3)
#lz[["zs.stat"]] <- s_rescale(lz[["z.final"]], factorv = meanv)
#lz[["zs.rand"]] <- s_rescale(lz[["z.final"]], meanv = meanv, sdv = sdv)
sef.fname <- "sef-markers_ct-treg_stransform-expt_dlpfc-ro1.rda"
sef.fpath <- file.path(proj.dpath, 
                       "outputs/03_test-stransform", 
                       sef.fname)
sef <- get(load(file.path(here(), sef.fpath)))
#lexpt <- lz[c("z.final", "zs.stat", "zs.rand")]