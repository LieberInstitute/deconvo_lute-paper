#!/usr/bin/env R

# Author: Sean Maden, Stephanie Hicks
#
# RNAscope quality on shuffle experiment status.
#

# source("00_param.R")

#------
# load
#------

load("./data/09_quality/halo_all.Rdata")

load("./env/03_shuffle/00_fig3ab_script.RData")

#-----------
# helper functions
#-----------
get_halo_df_conf <- function(halo.table, data.type = "halo_all",
                             sample.id.colname = "Sample"){
  # get_halo_df_conf
  #
  # Author: Sean Maden
  # 
  # halo.table : a tibble.
  # data.type : user-specified data type.
  #
  #
  #
  #
  
  if(!is(halo.table, "tbl")){stop("Error: object needs to be of type tbl.")}
  if(!is(halo.table, "tbl_df")){stop("Error: object needs to be of type tbl.")}
  if(!is(halo.table, "data.frame")){stop("Error: object needs to be of type tbl.")}
  
  sample.id.vector <- as.character(unique(as.matrix(halo.table[,sample.id.colname])))
  
  if(data.type == "halo_all"){
    combo.label.vector <- c("Circle", "Star")
    combo.label.varname <- "Combo"
    do.call(cind, lapply(combo.label.vector, function(combo.label){
      
    }))
    
    
  } else{
    stop("error, unidentified data type.")
  }
  df.returns <- df.conf
  return(df.returns)
}



#----------
# map data
#----------
# 

df.conf <- get_halo_df_conf()

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


