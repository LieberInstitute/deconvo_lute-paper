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

conf_by_combo <- function(halo.processed.table,
                          combo.label.varname){
  # conf_by_combo
  #
  # get confidence by slide combination.
  #
  #
  #
  #
  #
  #
  
  
  combo.label.vector <- as.vector(
    unique(
      halo.processed.table[,combo.label.varname]))[[1]]
  
  # get conf by combo
  df.conf.bycombo <- as.data.frame(
    
    table(
      
      halo.processed.table$Confidence, 
      halo.processed.table$Combo, 
      halo.processed.table$Sample))
  
  
  df.conf.bycombod <- df.conf.bycombo[!df.conf.bycombo[,4]==0, seq(3)]
  colnames(df.conf.bycombod) <- c("confidence", "combo", "sample.id")
  
  lr <- list(df.conf.bycombod = df.conf.bycombod)
  
  return(lr)
  
  
}



conf_from_halodataobject <- function(
    halo.table, data.type = "halo_all",
                             
    sample.id.colname = "Sample",
                             
    combo.label.varname = "Combo",
                             
    combo.label.vector = c("Circle", "Star")
    
    ){
  # conf_from_halodataobject
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
  
  halo.processed.table <- halo.table
  
  if(data.type == "halo_all"){
    
    
    df.conf <- conf_by_combo(halo.processed.table,
                             combo.label.varname)$df.conf.bycombod
    
    
  } else{
    stop("error, unidentified data type.")
  }
  
  lr <- list(df.conf.bycombo = df.conf)
  
  return(lr)
}


conf_frequencies <- function(df.conf,
                             high.label = "High",
                             low.label = "Low",
                             middle.label = "OK"){
  # get combo quality frequency
  #
  #
  # df.conf : Confidence table with standard colnames (confidence, combo, sample.id).
  # 
  # status = "both low" # consensus low/below quality thresh
  # status = "both high" # consensus high/exceeds quality thresh
  # status = "both OK" # consensus ok/meets quality thresh
  # status = "missmatch.lowhigh" # 
  # status = "missmatch.okhigh" # cross slide variances?
  # stats = "missmatch.oklow" # moderate/marginal
  #
  # is.high <- df.conf$confidence==high.label
  # is.low <- df.conf$confidence==low.label
  # is.middle <- df.conf$confidence==middle.label
  #
  # is.high.and.middle <- is.high & is.middle
  # is.high.and.low <- is.high & is.low
  # is.low.and.middle <- is.low & is.middle
  #
  # note:
  # we may retain the at or above quality samples
  # in other words, removes low/high, removes oklow, removes both low
  #
  # RETURNS:
  # TABLE OF 
  # (ROWS) UNIQUE SAMPLE IDS
  # (COLUMNS) STATUS TYPE (ENTRIES: T/F)
  #
  #
  #
  
  df.conf.all <- df.conf
  sample.id.vector <- unique(df.conf.slide$sample.id)
  colnames.iter <- c("sample.id", "is.high", "is.low", "is.middle",
                     "is.high.and.middle", "is.high.and.low", "is.low.and.middle")
  
  list.combo.return <- lapply(
    unique(df.conf$combo), 
    function(combo){
      df.conf <- df.conf[df.conf$combo==combo,]
      is.high <- df.conf$confidence==high.label
      is.low <- df.conf$confidence==low.label
      is.middle <- df.conf$confidence==middle.label
      is.high.and.middle <- is.high & is.middle
      is.high.and.low <- is.high & is.low
      is.low.and.middle <- is.low & is.middle
      list.return <- lapply(sample.id.vector, 
                            function(sample.id){
        dff.conf <- df.conf[df.conf$sample.id==sample.id,]
        is.high <- dff.conf$confidence==high.label
        is.low <- dff.conf$confidence==low.label
        is.middle <- dff.conf$confidence==middle.label
        is.high.and.middle <- is.high & is.middle
        is.high.and.low <- is.high & is.low
        is.low.and.middle <- is.low & is.middle
        return(
          c(as.character(sample.id),
             is.high,
             is.low,
             is.middle,
             is.high.and.middle,
             is.high.and.low,
             is.low.and.middle)
          )
      })
      df.return <- 
        as.data.frame(
        do.call(rbind, 
                lapply(
                  list.return, 
                  function(item){item})))
      colnames(df.return) <- colnames.iter
      df.return$combo <- combo
      return(df.return)
    }
    )
  
  # get wide df
  df.wide <- as.data.frame(do.call(cbind, 
                                   lapply(list.combo.return, function(item){
                                     combo <- unique(item$combo)
                                     colnames(item) <- paste0(colnames(item), ".", combo)
                                     item
                                    })))
  
  df.wide$is.high.consensus <- 
    as.logical(df.wide$is.high.Circle) & as.logical(df.wide$is.high.Star)
  df.wide$is.low.consensus <- 
    as.logical(df.wide$is.low.Circle) & as.logical(df.wide$is.low.Star)
  df.wide$is.middle.consensus <- 
    as.logical(df.wide$is.middle.Circle) & as.logical(df.wide$is.middle.Star)
  
  df.wide$is.high.and.low <- 
    as.logical(df.wide$is.high.Circle) & as.logical(df.wide$is.low.Star) |
    as.logical(df.wide$is.low.Circle) & as.logical(df.wide$is.high.Star)
  
  df.wide$is.high.and.middle <- 
    as.logical(df.wide$is.high.Circle) & as.logical(df.wide$is.middle.Star) |
    as.logical(df.wide$is.middle.Circle) & as.logical(df.wide$is.high.Star)
  
  df.wide$is.low.and.middle <- 
    as.logical(df.wide$is.high.Circle) & as.logical(df.wide$is.low.Star) |
    as.logical(df.wide$is.low.Circle) & as.logical(df.wide$is.high.Star)
  
  df.combo.wide <- df.wide
  
  lr <- list(list.combo.return = list.combo.return,
             df.combo.wide = df.combo.wide)
  
  return(lr)
  
  
}

