cellCategorySummaries <- function(cellTypeVector, typeLabelSubset=5){
  totalCell <- length(cellTypeVector)
  
  isTcellString <- "^T .*|^T\\..*"
  isBcellString <- "^B .*|^B\\..*"
  
  isTcell <- grepl(isTcellString, cellTypeVector)
  isBcell <- grepl(isBcellString, cellTypeVector)
  isOther <- !(isTcell|isBcell)
  
  numTcell <- length(which(isTcell))
  numBcell <- length(which(isBcell))
  numOther <- length(which(isOther))
  
  fractionTcell <- round(numTcell/totalCell, 3)
  fractionBcell <- round(numBcell/totalCell, 3)
  fractionOther <- round(numOther/totalCell, 3)
  
  tCellLabels <- paste0(cellTypeVector[isTcell],collapse=";")
  bCellLabels <- paste0(cellTypeVector[isBcell],collapse=";")
  otherCellLabels <- paste0(cellTypeVector[isOther],collapse=";")
  
  if(typeLabelSubset){
    tCellLabels <- tCellLabels[seq(typeLabelSubset)]
    bCellLabels <- bCellLabels[seq(typeLabelSubset)]
    otherCellLabels <- otherCellLabels[seq(typeLabelSubset)]
  } else{}
  
  dfSummary <- data.frame(
    type = c("T-cell", "B-cell", "Other (not T- or B-cell)"),
    typeCount = c(numTcell, numBcell, numOther),
    typeFraction = c(fractionTcell, fractionBcell, fractionOther),
    typeLabels = c(tCellLabels, bCellLabels, otherCellLabels)
  )
  rownames(dfSummary)
  
  return(dfSummary)
}
