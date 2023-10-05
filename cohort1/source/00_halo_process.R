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
  #
  #
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
  unique.combo.vector <- unique(df.conf$combo)
  colnames.iter <- c("sample.id", "is.high", "is.low", "is.middle",
                     "is.high.and.middle", "is.high.and.low", "is.low.and.middle")
  
  list.combo.return <- lapply(
    unique.combo.vector, 
    function(combo){
      
      df.conf <- df.conf.all[df.conf.all$combo==combo,]
      sample.id.vector.iter <- df.conf$sample.id # all samples
      
      is.high <- df.conf$confidence==high.label
      is.low <- df.conf$confidence==low.label
      is.middle <- df.conf$confidence==middle.label
      is.high.and.middle <- is.high & is.middle
      is.high.and.low <- is.high & is.low
      is.low.and.middle <- is.low & is.middle
      
      item.na <- rep("NA", 7) # make null set, define NA terms
      
      list.combo.return <- lapply(sample.id.vector, 
                            function(sample.id){
                              
        message(sample.id)
        item.na.iter <- item.na
        item.na.iter[1] <- as.character(sample.id)
        
        # if slide available for sample, parse options, otherwise pass NA
        if(sample.id %in% sample.id.vector.iter){
          df.conf.iter <- df.conf[df.conf$sample.id==sample.id,]
          
          is.high <- df.conf.iter$confidence==high.label
          is.low <- df.conf.iter$confidence==low.label
          is.middle <- df.conf.iter$confidence==middle.label
          
          is.high.and.middle <- is.high & is.middle
          is.high.and.low <- is.high & is.low
          is.low.and.middle <- is.low & is.middle
          
          item.na.iter[2] <- is.high
          item.na.iter[3] <- is.low
          item.na.iter[4] <- is.middle
          item.na.iter[5] <- is.high.and.middle
          item.na.iter[6] <- is.high.and.low
          item.na.iter[7] <- is.low.and.middle
        }
        
        return(item.na.iter)
        
      })
      
      df.return <- 
        as.data.frame(
          do.call(rbind, 
                  lapply(
                    list.combo.return, 
                    function(item){item})))
      
      colnames(df.return) <- colnames.iter
      df.return$combo <- combo
      
      return(df.return)
    })
  
  # get wide df
  matrix.wide <- do.call(
    cbind, lapply(list.combo.return, function(item){
                           # message(item)
                           combo <- unique(item$combo)
                           colnames(item) <- 
                             paste0(colnames(item), ".", combo)
                           item
                          }))
  df.wide <- as.data.frame(matrix.wide)
  
  # get consensus/harmonized annotations ACROSS slides
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
  
  # has slides condition
  condition.has.all.slides <- 
    df.wide$is.high.Circle=="NA" | df.wide$is.high.Star=="NA"
  df.wide$has.all.slides <- !condition.has.all.slides
  
  # return
  lr <- list(
    list.combo.return = list.combo.return,
    df.combo.wide = df.wide
  )
  return(lr)
  
  
}

