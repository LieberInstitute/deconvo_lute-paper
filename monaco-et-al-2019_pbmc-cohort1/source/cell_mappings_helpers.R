#!/usr/bin/env R

# Author: Sean Maden
#
#
#
#
#
#

cellLabelMappings <- function(
    vectorCellTypeMap, vectorCellTypeStart = NULL, dfToMap = NULL, 
    returnType = c("list", "df"), summaryOperation = c("mean", "median", "sum")){
  #
  #
  #
  #
  # uses regex to map labels
  # returnType : whether to return list or single df (mapped values)
  # summaryOperation : how values were summarized
  #
  #
  
  # get mappings
  if(is(vectorCellTypeStart, "NULL")){
    mappingsTable <- data.frame(celltype1 = colnames(dfToMap))
  } else{
    mappingsTable <- data.frame(celltype1 = vectorCellTypeStart)
  }
  
  mappingsTable$celltype2 <- "NA"
  for(type in vectorCellTypeMap){
    message(type)
    filterCelltype1 <- 
      which(grepl(paste0("^", type, ".*"), mappingsTable$celltype1))
    if(length(filterCelltype1)>0){
      mappingsTable[filterCelltype1,]$celltype2 <- type
    }
  }
  
  # map terms to reference
  dfMapped <- "NULL"
  if(is(dfToMap, "NULL")){
    message("Skipping mapping summaries")
  } else{
    dfMapped <- dfToMap
    colnames(dfMapped) <- unlist(lapply(colnames(dfMapped),function(cname){
      mappingsTable[mappingsTable[,1]==cname,2]
    }))
    vectorNewTypes <- unique(mappingsTable[,2])
    dfMapped <- do.call(cbind, lapply(vectorNewTypes, function(type){
      filterType <- which(
        colnames(dfToMap) %in% mappingsTable[mappingsTable[,2]==type,1])
      if(length(filterType)>1){
        if(summaryOperation=="mean"){
          rowMeans(dfToMap[,filterType])
        } else if(summaryOperation=="median"){
          rowMedians(dfToMap[,filterType])
        } else{
          rowSums(dfToMap[,filterType])
        }
      } else{
        dfToMap[,filterType]
      }
    }))
    colnames(dfMapped) <- vectorNewTypes
  }
  
  # parse return options
  if(returnType == "list"){
    returnObject <- list(
      mappingsTable = mappingsTable,
      dfMapped = dfMapped
    )
  } else if(returnType == "df" & is(dfMapped, "NULL")){
    returnObject <- mappingsTable
  } else{
    returnObject <- dfMapped
  }
  return(returnObject)
}