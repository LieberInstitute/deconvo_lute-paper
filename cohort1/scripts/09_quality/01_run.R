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

list.dfp.tall <- list(
  
  dfp.tall = dfp.tall,
  dfp.tall.high = dfp.tall.high,
  dfp.tall.low = dfp.tall.low
  
)

df.conf.wide <- df.conf.freq$df.combo.wide

list.dfp.tall.append.quality <- lapply(list.dfp.tall, function(item){
  dfp.tall.iter <- item
  
  dfp.tall.iter$has.all.slides <- 
    dfp.tall.iter$is.high.and.low <- 
    dfp.tall.iter$is.high.and.middle <-
    dfp.tall.iter$is.low.and.middle <-
    "NA"
  dfp.tall.iter$is.high.consensus <- 
    dfp.tall.iter$is.low.consensus <- 
    dfp.tall.iter$is.middle.consensus <-
    "NA"
  
  for(sample.id in df.conf.wide$sample.id.Circle){
    
    
    message(sample.id)
    df.conf.wide.iter <- 
      df.conf.wide[df.conf.wide$sample.id.Circle==sample.id,,drop = F]
    
    
    if(sample.id %in% dfp.tall.iter$sample.id){
      iter.filter.condition <- dfp.tall.iter$sample.id==sample.id
      
      dfp.tall.iter[iter.filter.condition,]$has.all.slides <- 
        df.conf.wide.iter$has.all.slides
      dfp.tall.iter[iter.filter.condition,]$is.high.and.low <- 
        df.conf.wide.iter$is.high.and.low
      dfp.tall.iter[iter.filter.condition,]$is.high.and.middle <- 
        df.conf.wide.iter$is.high.and.middle
      dfp.tall.iter[iter.filter.condition,]$is.low.and.middle <- 
        df.conf.wide.iter$is.low.and.middle
      dfp.tall.iter[iter.filter.condition,]$is.high.consensus <- 
        df.conf.wide.iter$is.high.consensus
      dfp.tall.iter[iter.filter.condition,]$is.low.consensus <- 
        df.conf.wide.iter$is.low.consensus
      dfp.tall.iter[iter.filter.condition,]$is.middle.consensus <- 
        df.conf.wide.iter$is.middle.consensus
    }
    
  }
  
  return(dfp.tall.iter)
  
})

names(list.dfp.tall.append.quality) <- names(list.dfp.tall)

#---------
# save
#---------

# save image
save.image("./env/09_quality/01_run_script.RData")


