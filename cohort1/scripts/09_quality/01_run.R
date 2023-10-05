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

df.conf.wide <- df.conf.freq$df.combo.wide
dfp.tall$has.all.slides <- 
  dfp.tall$is.high.and.low <- 
  dfp.tall$is.high.and.middle <-
  dfp.tall$is.low.and.middle <-
  "NA"
dfp.tall$is.high.consensus <- 
  dfp.tall$is.low.consensus <- 
  dfp.tall$is.middle.consensus <-
  "NA"

for(sample.id in df.conf.wide$sample.id.Circle){
  df.conf.wide.iter <- df.conf.wide[df.conf.wide$sample.id.Circle==sample.id,]
  dfp.tall$has.all.slides <- df.conf.wide.iter$has.all.slides
  dfp.tall$is.high.and.low <- df.conf.wide.iter$is.high.and.low
  dfp.tall$is.high.and.middle <- df.conf.wide.iter$is.high.and.middle
  dfp.tall$is.low.and.middle <- df.conf.wide.iter$is.low.and.middle
  dfp.tall$is.high.consensus <- df.conf.wide.iter$is.high.consensus
  dfp.tall$is.low.consensus <- df.conf.wide.iter$is.low.consensus
  dfp.tall$is.middle.consensus <- df.conf.wide.iter$is.middle.consensus
}

#---------
# save
#---------

# save image
save.image("./env/09_quality/01_run_script.RData")


