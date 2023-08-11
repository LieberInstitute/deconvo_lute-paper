rnascope_cell_info <- function(df.rn, sample.id = NULL, k.type = NULL, cell.types = NULL){
  #
  #
  # get cell info data from summary df
  #
  # example:
  # df.rn <- mae[["df.cellstat.rnascope"]]
  # ct.info <- rnascope_cell_info(df.rn, "k2", c("glial", "neuron"))
  #
  filter.vector <- colnames(df.rn)
  filter.condition.vector <- rep(TRUE, length(filter.vector))
  ct.string <- paste0(cell.types, collapse = "|")
  if(!is(ct.string, "NULL")){
    filter.condition.vector <- filter.condition.vector & grepl(ct.string, filter.vector)
  }
  if(!is(sample.id, "NULL")){
    filter.condition.vector <- filter.condition.vector & grepl(sample.id, filter.vector)
  }
  if(!is(k.type, "NULL")){
    filter.condition.vector <- filter.condition.vector & grepl(k.type, filter.vector)
  }
  df.ct <- df.rn[,filter.condition.vector]
  return(df.ct)
}

df.true.list <- function(df.rn, sample.id.vector, k.type, cell.types, info = "true_proportion"){
  #
  #
  # example:
  # sample.id.vector <- unique(y.unadj$batch.id2)
  # list.dftrue <- df.true.list(df.rn, sample.id.vector, "k2", c("glial", "neuron"))
  # names(list.dftrue) <- sample.id.vector
  #
  #
  
  list.dfinfo <- lapply(sample.id.vector, function(sample.id){
    rnascope_cell_info(df.rn, sample.id, k.type, cell.types)
  })
  list.dfinfo <- lapply(list.dfinfo, function(dfinfo){
    as.data.frame(dfinfo["true_proportion",,drop=F])
  })
  names(list.dfinfo) <- sample.id.vector
  return(list.dfinfo)
}
