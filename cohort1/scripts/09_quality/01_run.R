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

# df.conf <- get_halo_df_conf()

# 
as.data.frame(table(halo_all$Confidence, halo_all$SAMPLE_ID))
df.conf <- as.data.frame(table(halo_all$Confidence, halo_all$Sample))
dfp.tall$confidence <- NA
for(sample.id in df.qual[,2]){
  qual.iter <- df.conf[df.conf[,2]==sample.id,1]
  print(qual.iter)
  dfp.tall[dfp.tall$sample.id==sample.id,]$confidence <- 
    as.character(unique(df.conf[df.conf[,2]==sample.id, 1]))
}


#---------
# save
#---------

# save image
save.image("./env/09_quality/01_run_script.RData")


