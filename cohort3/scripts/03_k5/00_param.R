#------------------
# helper functions
#------------------

# FROM EXPRESSION, GET SIMPLIFIED THRESHOLD PLOT
getQuantileTablesFromReferenceExpression <- function(
    expressionTable, heatmapNumericScaleName, quantileThresholdIndex=4, 
    quantileSeq=seq(0,1,0.25)
){
  
  # boolean matrix
  booleanTable <- apply(tpmReference, 2, function(ci){
    quantileThreshold <- 
      quantile(ci, probs=quantileSeq)[quantileThresholdIndex]
    labelColumnValues <- ifelse(
      ci >= quantileThreshold, TRUE, FALSE)
    return(labelColumnValues)
  })
  # numeric matrix
  numericTable <- matrix(
    as.numeric(booleanTable), nrow = nrow(booleanTable))
  # format
  colnames(numericTable) <- 
    colnames(booleanTable) <- colnames(expressionTable)
  rownames(numericTable) <- 
    rownames(booleanTable) <- rownames(expressionTable)
  
  # heatmaps
  heatmapNumeric <- Heatmap(
    numericTable, name = heatmapNumericScaleName)
  
  # return
  returnList <- list(
    booleanTable = booleanTable, 
    numericTable = numericTable,
    heatmapNumeric = heatmapNumeric
  )
  return(returnList)
}

# FROM EXPRESSION GET PCA

# 