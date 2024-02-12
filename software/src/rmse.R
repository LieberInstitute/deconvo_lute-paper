#!/usr/bin/env R

# Author: Sean Maden
#
# Get RMSE calculations compatiable with test datasets.
#
#
#

labelsVectorFromFilter <- function(filterLabels, filterConditions,
                                   labels = c("sample.id", "condition.id", "replicate.id")){
  # labelsVectorFromFilter
  #
  # Bind filter term labels from filter list. Used to obtain condition labels from filter list.
  # See also applyRmse().
  #
  
  item1 <- filterLabels
  item1Label <- labels[item1]
  
  item2 <- filterConditions
  item2Unlist <- unlist(strsplit(item2, ";"))
  
  labels <- paste0(rep(item1Label, length(item2)), ":", item2Unlist)
  
  labelsVector <- c()
  incrementBind <- length(item2Unlist)/length(item2)
  start <- 1
  for(i in seq(item2)){
    end <- start + incrementBind - 1
    labelsVector <- c(labelsVector,
                      paste(labels[start:end], collapse = ";"))
    start <- end + 1
  }
  labelsVector
}

labelsVectorFromFilter <- function(filterLabels, filterConditions,
                                   labels = c("sample.id", "condition.id", "replicate.id")){
  # labelsVectorFromFilter
  #
  # bind filter term labels from filter list
  #
  
  item1 <- filterLabels
  item1Label <- labels[item1]
  
  item2 <- filterConditions
  item2Unlist <- unlist(strsplit(item2, ";"))
  
  labels <- paste0(rep(item1Label, length(item2)), ":", item2Unlist)
  
  labelsVector <- c()
  incrementBind <- length(item2Unlist)/length(item2)
  start <- 1
  for(i in seq(item2)){
    end <- start + incrementBind - 1
    labelsVector <- c(labelsVector,
                      paste(labels[start:end], collapse = ";"))
    start <- end + 1
  }
  labelsVector
}

dataSubsetOnFilter <- function(data, subset, filter){
  # dataSubsetOnFilter
  # data, valid rmse data table
  # subset, subset mapping
  # filter, filter vector
  #
  colNames <- c("sample.id", "condition.id", "replicate.id")
  fv <- unlist(strsplit(as.character(filter), ";"))
  conditionVector <- rep(TRUE, nrow(data))
  filterIndex <- 1
  for(s in subset){
    conditionVector <- 
        conditionVector & data[,colNames[s]]==fv[filterIndex]
    filterIndex <- filterIndex + 1
  }
  data[conditionVector,]
}

listFilterGroups <- function(rmseData, 
                             colNamesVector = 
                               c("sample.id", "condition.id", "replicate.id"),
                             sepSymbol = ";"){
	# listFilterGroups
	# 
	#
  rmseDataSubset <- rmseData[,colNamesVector]
  whichSubgroup <- which(colnames(rmseDataSubset)=="sample.id")
  whichCondition <- which(colnames(rmseDataSubset)=="condition.id")
  whichReplicate <- which(colnames(rmseDataSubset)=="replicate.id")
  
  listFilterData <- list()
  subsetVectorList <- list(
    c(whichSubgroup),
    c(whichCondition),
    c(whichReplicate),
    c(whichReplicate, whichCondition)
  )
  for(subsetIter in subsetVectorList){
    if(length(subsetIter)==1){
      listFilterData[[length(listFilterData) + 1]] <-
      list(
        subset = subsetIter,
        filter = paste0(unique(rmseDataSubset[,subsetIter]))
      )
    } else{
      listFilterData[[length(listFilterData) + 1]] <-
      list(
        subset = subsetIter,
        filter = paste0(unique(rmseDataSubset[,subsetIter[1]]),
                        sepSymbol,
                        rep(
                          unique(rmseDataSubset[,subsetIter[2]]),
                          each = length(
                            unique(rmseDataSubset[,subsetIter[1]])
                          )
                        ))
      )
    }
  }
  return(listFilterData)
}

returnFilteredProportionsData <- function(
    proportionsData, listFilterData){
  # returnFilteredProportionsData
  #
  # Get list of filtered datasets.
  listSubsetAll <- lapply(listFilterData, function(filterParam){
  listDataSubset <- lapply(filterParam[[2]], function(f){
    dataSubset <- dataSubsetOnFilter(proportionsData,
                       subset = filterParam[[1]],
                       filter = f)
  })
  return(listDataSubset)
  })
  listSubsetUnnest <- unlist(listSubsetAll, recursive = FALSE)
}

rmseTallTableFromResult <- function(rmseResult){
  # rmseTallTableFromResult
  #
  # Gets tall table of RMSE results.
  #
  lfp <- rmseResult$listFiltersParams
conditionsVector <- lapply(lfp, function(item){
  labelsVectorFromFilter(item$subset, item$filter)
}) |> unlist()
rmseResultVector <- rmseResult$returnRmse$values |> unlist()
rmseTableTall <- data.frame(condition = conditionsVector,
                            rmseResult = rmseResultVector)
  return(
    rmseTableTall
  )
}

absErrorsFromTypes <- function(proportion1, proportion2){
  # absErrorsFromTypes
  #
  # takes absolute of 1:1 difference from arguments
  #
  # use to check errors
  # check errors
  #errorsVectorCheck1 <- filterDataTypes[,level1Filter][1]-
  #  filterDataTypes[,level2Filter][1]
  #errorsVectorCheck2 <- absErrorsFromTypes(
  #    filterDataTypes[,level1Filter],
  #    filterDataTypes[,level2Filter])[1,1]
  #errorsVectorCheck1[1]==errorsVectorCheck2[1,1]
  #
  #
  abs(proportion1-proportion2)
}

filterDataRmseList <- function(filterDataList, cellTypesVector){
  lapply(filterDataList, function(filterData){
    filterData$rmse <- 
      rmseParam_rootMeanSquaredError(filterData, cellTypesVector)
    return(filterData)
  })
}

formatRmseReturn <- function(filterDataList, cellTypesVector = NULL){
  valuesByFilter <- rmseValuesByFilter(filterDataList, cellTypesVector)
  valuesDataList <- filterDataRmseList(filterDataList, cellTypesVector)
  tableValues <- 
    do.call(rbind, valuesDataList) |> as.data.frame()
  tableValuesSparse <- valuesByFilter |> 
    unlist() |> as.data.frame()
    
  returnList <- list(
    values = valuesByFilter,
    list = valuesDataList,
    tableValues = tableValues,
    tableValuesSparse
  )
}

rmseParam_rootMeanSquaredError <- function(
    filterData, 
    cellTypesVector = c("celltype1", "celltype2"),
    levelsVector = c("known", "predicted"),
    detail = FALSE){
  # rmseParam_rootMeanSquaredError
  #
  # filterData, filtered valid rmse dataset
  # cellTypesVector, vector of celltypes in colnames start
  # levelsVector, vector of levels in colnames end
  #
  # returns vector of rmse values of dim CELLTYPES
  # if detail true, return env components as list.
  #
  #
  #
  #
  #
  
  typesFilter <- grepl(
    paste0("^",cellTypesVector, ".*", collapse = "|"),
    colnames(filterData))
  typesFilter <- typesFilter &
    grepl(
      paste0(".*",levelsVector, "$", collapse = "|"),
      colnames(filterData))
  filterDataTypes <- filterData[,typesFilter]
  
  level1Filter <- 
    grepl(paste0(".*", levelsVector[1], "$"),
          colnames(filterDataTypes))
  level2Filter <- 
    grepl(paste0(".*", levelsVector[2], "$"),
          colnames(filterDataTypes))
  
  # RMSE OF TYPES BY GROUP * CELLTYPE
  errorsVector <- absErrorsFromTypes(
    filterDataTypes[,level1Filter],
    filterDataTypes[,level2Filter]) |> unlist()
  errorsSq <- errorsVector^2
  denomMeanOperator <- length(errorsSq) # length(unique(cellTypesVector))
  # get sum of squared errors vector
  numGroups <- nrow(filterDataTypes[,level1Filter,drop=F])
  sumErrorsByGroup <- c()
  seqIndexErrors <- 
    seq(from = 1, to = length(errorsSq), 
        by = length(errorsSq)/denomMeanOperator)
  for(type in cellTypesVector){
    startIndex <- seqIndexErrors[which(cellTypesVector==type)]
    stopIndex <- startIndex + numGroups - 1
    newSum <- sum(errorsSq[startIndex:stopIndex])
    sumErrorsByGroup <- c(sumErrorsByGroup, newSum)
  }
  # test
  # sum(errorsSq[1:3])==sumErrorsByGroup[1]
  # sum(errorsSq[4:6])==sumErrorsByGroup[2]
  errorsSqMean <- sum(errorsSq)/denomMeanOperator # mean(sumErrorsByGroup)
  rmseReturn <- sqrt(errorsSqMean)
  if(detail){
    rmseReturnList <- list(
      denomMeanOperator = denomMeanOperator,
      sumErrorsByGroup = sumErrorsByGroup,
      errorsSqMean = errorsSqMean,
      rmseReturn = rmseReturn
    )
    return(rmseReturnList)
  } else{
    return(rmseReturn)
  }
  return(FALSE)
}

rmseValuesByFilter <- function(filterDataList, cellTypesVector = NULL){
  lapply(filterDataList, function(item){
    if(nrow(item)>0){
      rmseParam_rootMeanSquaredError(item, cellTypesVector)
      }
    })
}

applyRmse <- function(rmseData, listFilterData = NULL, cellTypesVector = NULL, detail = FALSE){
  # applyRmse
  #
  # rmseData
  # listFilterData
  # cellTypesVector
  # detail
  #
  if(is(listFilterData, "NULL")){
    listFiltersParams <- listFilterGroups(rmseData)
  } else{
    listFiltersParams <- listFilterData
  }
  filterDataList <- 
    returnFilteredProportionsData(rmseData, listFiltersParams)
  rmseOut <- formatRmseReturn(filterDataList, cellTypesVector)
  
  returnList <- list(
    "filterDataList" = filterDataList, 
    "listFiltersParams" = listFiltersParams, 
    "returnRmse" = rmseOut
  )
  returnList[["returnRmseTall"]] <- rmseTallTableFromResult(returnList)
  return(returnList)
}