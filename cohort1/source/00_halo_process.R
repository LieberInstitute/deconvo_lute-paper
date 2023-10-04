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



get_halo_df_conf <- function(halo.table, data.type = "halo_all",
                             sample.id.colname = "Sample",
                             combo.label.varname = "Combo"){
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
  
  combo.label.vector <- c("Circle", "Star")
  combo.label.varname <- "Combo"
  
  if(data.type == "halo_all"){
    
    
    # get conf by combo
    df.conf.bycombo <- as.data.frame(
      table(
        halo.table$Confidence, halo.table$Combo, halo.table$Sample))
    df.conf.bycombod <- dft[!dft[,4]==0, seq(3)]
    colnames(df.conf.bycombod) <- c("confidence", "combo", "sample.id")
    
    
  } else{
    stop("error, unidentified data type.")
  }
  
  # get combo quality frequency
  # status = "both low" # consensus low/below quality thresh
  # status = "both high" # consensus high/exceeds quality thresh
  # status = "both OK" # consensus ok/meets quality thresh
  # status = "missmatch.lowhigh" # 
  # status = "missmatch.okhigh" # cross slide variances?
  # stats = "missmatch.oklow" # moderate/marginal
  #
  # note:
  # we may retain the at or above quality samples
  # in other words, removes low/high, removes oklow, removes both low
  
  
  lr <- list(df.conf.bycombo = df.conf.bycombod)
  
  df.returns <- df.conf
  
  return(df.returns)
}