#!/usr/bin/env R

# Author: Sean Maden, Stephanie Hicks
#
# RNAscope quality on shuffle experiment status.
#
#
#
#
#
#

source("./source/00_halo_process.R")

#------
# load
#------



load("./data/09_quality/halo_all.Rdata")

load("./env/03_shuffle/00_fig3ab_script.RData")


#----------
# map data
#----------
# 

df.conf <- conf_from_halodataobject(halo_all)[[1]]

df.conf.freq <- conf_frequencies(df.conf)

map_halo_conf_df <- function(df, df.conf, sample.id.variable = "sample.id"){
  
}


#---------
# save
#---------

# save image
save.image("./env/09_quality/01_run_script.RData")


