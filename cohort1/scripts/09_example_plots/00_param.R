#!/usr/bin/env R

# Author: Sean Maden
#
# Functions for example plot generation.
#
#
#
#

libv <- c("ggplot2", "reshape2", "gridExtra", "lute")
sapply(libv, library, character.only = T)

#-------------------------

# Single example functions
 
#-------------------------

singleValueTestVariables <- function(cellScaleFactorsStart = 2,
                                     cellScaleFactorOffTypeValue = 10,
                                     trueProportionValue = 0.8,
                                     markerExpressionStart = 3,
                                     cellScaleFactorNew = 10,
                                     bulkExpressionValue = 0.8){
  # singleValueTestVariables
  #
  # Run to get example valuesList object.
  #
  #
  #
  #
  
  valuesList <- list(
    cellScaleFactorsStart = cellScaleFactorsStart,
    cellScaleFactorOffTypeValue = cellScaleFactorOffTypeValue,
    trueProportionValue = trueProportionValue,
    markerExpressionStart = markerExpressionStart,
    cellScaleFactorNew = cellScaleFactorNew,
    bulkExpressionValue = bulkExpressionValue
  )
  
  valuesList <- parseExampleStartValues(valuesList)
  
  return(valuesList)
}

parseExampleStartValues <- function(valuesList){
  # parseExampleStartValues
  #
  # Get predictions for starting values, and append to valuesList
  #
  #
  #
  #
  
  # get scale factor sets
  cellScaleFactorsStart <- c(valuesList[["cellScaleFactorsStart"]], 
                             valuesList[["cellScaleFactorOffTypeValue"]])
  cellScaleFactorsNull <- c(1, 1)
  names(cellScaleFactorsStart) <- 
    names(cellScaleFactorsNull) <- c("type1", "type2")
  
  # get expression matrices
  matrixValues <- c(valuesList[["markerExpressionStart"]], 0, 0, 
                    valuesList[["markerExpressionStart"]])
  zrefExample <- matrix(matrixValues, nrow = 2)
  colnames(zrefExample) <- c("type1", "type2")
  bulkExpressionExample <- matrix(
    rep(valuesList[["bulkExpressionValue"]], 2), ncol = 1)
  rownames(zrefExample) <- rownames(bulkExpressionExample) <- 
    paste0("gene", seq(nrow(zrefExample)))
  
  newParamStart <- 
    nnlsParam(bulkExpressionExample, zrefExample, cellScaleFactorsStart) |>
    deconvolution()
  newParamNull <- 
    nnlsParam(bulkExpressionExample, zrefExample, cellScaleFactorsNull) |>
    deconvolution()
  predictedProportionsStart <- newParamStart@predictionsTable[[1]]
  predictedProportionsNull <- newParamNull@predictionsTable[[1]]
  biasStart <- valuesList[["trueProportionValue"]]-predictedProportionsStart
  biasNull <- valuesList[["trueProportionValue"]]-predictedProportionsNull
  errorStart <- abs(biasStart)
  errorNull <- abs(biasNull)
  
  # return
  valuesList[["deconvoResultsStart"]] <- newParamStart
  valuesList[["bulkExpressionExample"]] <- bulkExpressionExample
  valuesList[["zrefExample"]] <- zrefExample
  valuesList[["predictedProportionsStart"]] <- predictedProportionsStart
  valuesList[["predictedProportionsNull"]] <- predictedProportionsNull
  valuesList[["biasStart"]] <- biasStart
  valuesList[["biasNull"]] <- biasNull
  valuesList[["errorStart"]] <- errorStart
  valuesList[["errorNull"]] <- errorNull
  return(valuesList)
}

singleValueExample <- function(valuesList, conditionLabel = ""){
  # singleValueExample
  #
  # Get single point values example. 
  #
  # @param valuesList Output from singleValueTestVariables().
  #
  #
  #
  #
  
  
  
  changeNew <- valuesList[["cellScaleFactorNew"]]-
    valuesList[["cellScaleFactorStart"]]
  labelNew <- 
    paste0("cellScaleFactor = ",valuesList[["cellScaleFactorNew"]],")")
  
  # get simulation results
  cellScaleFactorsNew <- c(valuesList[["cellScaleFactorNew"]], 
                           valuesList[["cellScaleFactorOffTypeValue"]])
  names(cellScaleFactorsNew) <- c("type1", "type2")
  newParamNew <- 
    nnlsParam(
      valuesList[["bulkExpressionExample"]], valuesList[["zrefExample"]], 
      cellScaleFactorsNew) |>
    deconvolution()
  predictedProportionsNew <- newParamNew@predictionsTable[[1]]
  biasNew <- trueProportionValue-predictedProportionsNew
  errorNew <- abs(biasNew)
  
  # get changes
  cellScaleFactorChange <- valuesList[["cellScaleFactorNew"]]-
    valuesList[["cellScaleFactorsStart"]]
  proportionChange <- predictedProportionsNew-
    valuesList[["predictedProportionsStart"]]
  biasChange <- biasNew-valuesList[["biasStart"]]
  errorChange <- errorNew-valuesList[["errorStart"]]
  markerExpressionChange <- 
    (valuesList[["markerExpressionStart"]]/valuesList[["cellScaleFactorNew"]])-
    valuesList[["markerExpressionStart"]]
  
  # get plot data
  dfpNew <- data.frame(
    variable = 
      c("cellScaleFactor", "markerExpression", 
        "predictedProportion", "bias", "error"),
    value = 
      c(cellScaleFactorChange, markerExpressionChange, 
        proportionChange, biasChange, errorChange)
  )
  dfpNew$conditionLabel <- conditionLabel
  dfpNew$Change <- ifelse(dfpNew$value > 0, "Increase", "Decrease")
  dfpNew$variable <- factor(dfpNew$variable, 
                            levels = c("cellScaleFactor", "markerExpression",
                                       "predictedProportion", "bias", "error"))
  
  plot2 <- ggplot(dfpNew, aes(x = variable, y = value, fill = Change)) + 
    geom_bar(stat="identity", color = "black") + theme_bw() +
    ylab("Change (New - Old)") + facet_wrap(~conditionLabel, nrow = 1) + 
    geom_hline(yintercept = 0) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    ggtitle(plotTitleStringChanges) +
    scale_fill_manual(breaks = changeLevelsVector, 
                      values=c("dodgerblue", "gold"))
  
  returnList <- list(
    valuesList = valuesList,
    deconvoResult = newParamNew,
    plotData = dfpNew,
    plot = plot2
  )
  
  return(returnList)
  
}


valuesList <- singleValueTestVariables(
  cellScaleFactorsStart = 0, cellScaleFactorNew = 2, trueProportionValue = 0.1)
exampleResult <- singleValueExample(valuesList)
exampleResult$plot

#-----------------------------------

# Faceted/combined example functions

#-----------------------------------

source_test_values <- function(){
  cellScaleFactorStart = 2
  cellScaleFactorOffTypeValue = 10
  trueProportionValue = 0.8
  markerExpression = 3
  cellScaleFactorNewHigh = 3
  cellScaleFactorNewMidHigh = 2.5
  cellScaleFactorNewMidLow = 1.5
  cellScaleFactorNewLow = 1
  changeNull = 0
  yLabelStringChanges = "Change"
  plotTitleStringChanges = "Affect of scale change"
  
  
  valuesList <- list(
    cellScaleFactorStart = 2,
    cellScaleFactorOffTypeValue = 10,
    trueProportionValue = 0.8,
    markerExpression = 3,
    cellScaleFactorNewHigh = 3
  )
}

point_value_example <- function(cellScaleFactorStart = 2,
                                cellScaleFactorOffTypeValue = 10,
                                trueProportionValue = 0.8,
                                markerExpression = 3,
                                cellScaleFactorNewHigh = 3,
                                cellScaleFactorNewMidHigh = 2.5,
                                cellScaleFactorNewMidLow = 1.5,
                                cellScaleFactorNewLow = 1,
                                changeNull = 0,
                                yLabelStringChanges = "Change",
                                plotTitleStringChanges = "Affect of scale change"){
  # point_value_example
  #
  # Get example from simulation using provided point values
  #
  #
  #
  
  #-----------------------
  # parse simulation params
  #-----------------------
  
  changeHigh <- cellScaleFactorNewHigh-cellScaleFactorStart
  changeMidHigh <- cellScaleFactorNewMidHigh-cellScaleFactorStart
  changeMidLow <- cellScaleFactorNewMidLow-cellScaleFactorStart
  changeLow <- cellScaleFactorNewLow-cellScaleFactorStart
  labelHigh <- 
    paste0("Increase\n(cellScaleFactor = ",cellScaleFactorNewHigh,")")
  labelMidHigh <- 
    paste0("Moderate Inc.\n(cellScaleFactor = ",cellScaleFactorNewMidHigh,")")
  labelMidLow <- 
    paste0("Moderate Dec.\n(cellScaleFactor = ",cellScaleFactorNewMidLow,")")
  labelLow <- paste0("Decrease\n(cellScaleFactor = ",cellScaleFactorNewLow,")")
  
  changeLevelsVector <- c(labelHigh, labelMidHigh, labelMidLow, labelLow)
  
  #-----------------------
  # get simulation results
  #-----------------------
  # get expression matrices
  zrefExample <- matrix(c(markerExpression, 0, 0, markerExpression), nrow = 2)
  colnames(zrefExample) <- c("type1", "type2")
  bulkExpressionExample <- matrix(c(0.5,0.5), ncol = 1)
  
  # get scale factor sets
  cellScaleFactorsStart <- c(cellScaleFactorStart, cellScaleFactorOffTypeValue)
  cellScaleFactorsNull <- c(1, cellScaleFactorOffTypeValue)
  cellScaleFactorsHigh <- c(cellScaleFactorNewHigh, cellScaleFactorOffTypeValue)
  cellScaleFactorsMidHigh <- 
    c(cellScaleFactorNewMidHigh, cellScaleFactorOffTypeValue)
  cellScaleFactorsMidLow <- 
    c(cellScaleFactorNewMidLow, cellScaleFactorOffTypeValue)
  cellScaleFactorsLow <- c(cellScaleFactorNewLow, cellScaleFactorOffTypeValue)
  
  names(cellScaleFactorsStart) <- names(cellScaleFactorsNull) <-
    names(cellScaleFactorsHigh) <- names(cellScaleFactorsMidHigh) <- 
    names(cellScaleFactorsLow) <- names(cellScaleFactorsMidLow) <- 
    c("type1", "type2")
  
  newParam <- 
    nnlsParam(bulkExpressionExample, zrefExample, cellScaleFactorsStart) |>
    deconvolution()
  newParamNullScales <- 
    nnlsParam(bulkExpressionExample, zrefExample, cellScaleFactorsNull) |>
    deconvolution()
  newParamHigh <- 
    nnlsParam(bulkExpressionExample, zrefExample, cellScaleFactorsHigh) |>
    deconvolution()
  newParamMidHigh <- 
    nnlsParam(bulkExpressionExample, zrefExample, cellScaleFactorsMidHigh) |>
    deconvolution()
  newParamMidLow <- 
    nnlsParam(bulkExpressionExample, zrefExample, cellScaleFactorsMidLow) |>
    deconvolution()
  newParamLow <- 
    nnlsParam(bulkExpressionExample, zrefExample, cellScaleFactorsLow) |>
    deconvolution()
  
  predictedProportionsStart <- newParam@predictionsTable[[1]]
  predictedProportionsNull <- newParamNullScales@predictionsTable[[1]]
  predictedProportionsHigh <- newParamHigh@predictionsTable[[1]]
  predictedProportionsMidHigh <- newParamMidHigh@predictionsTable[[1]]
  predictedProportionsMidLow <- newParamMidLow@predictionsTable[[1]]
  predictedProportionsLow <- newParamLow@predictionsTable[[1]]
  
  biasStart <- trueProportionValue-predictedProportionsStart
  biasNull <- trueProportionValue-predictedProportionsNull
  biasHigh <- trueProportionValue-predictedProportionsHigh
  biasMidHigh <- trueProportionValue-predictedProportionsMidHigh
  biasMidLow <- trueProportionValue-predictedProportionsMidLow
  biasLow <- trueProportionValue-predictedProportionsLow
  errorStart <- abs(biasStart)
  errorNull <- abs(biasNull)
  errorHigh <- abs(biasHigh)
  errorMidHigh <- abs(biasMidHigh)
  errorMidLow <- abs(biasMidLow)
  errorLow <- abs(biasLow)
  
  #---------------------------------
  # make barplot with example values
  #---------------------------------
  dfp <- data.frame(
    cellScaleFactorStart = cellScaleFactorStart,
    markerExpression = markerExpression,
    trueProportion = trueProportionValue,
    predictedProportion = predictedProportionsStart,
    biasValue = biasStart,
    errorValue = errorStart
  )
  dfp <- melt(dfp)
  dfp$value <- round(dfp$value, 3)
  
  plot1 <- ggplot(dfp, aes(x = variable, y = value)) + 
    geom_bar(stat="identity", color = "black", fill = 'gray') + 
    theme_bw() + geom_hline(yintercept = 0) +
    theme(axis.text.x = element_text(angle=45,hjust=1),
          axis.title.x = element_blank()) +
    geom_text(aes(label = value), vjust = -1.5) +
    ylab("Value") + ylim(min(dfp$value-1), max(dfp$value+1))
  
  #-------------------
  # new direction plot
  #-------------------
  facetVector <- c(labelHigh, labelMidHigh, labelMidLow, labelLow)
  newScalesVector <- c(cellScaleFactorNewHigh, cellScaleFactorNewMidHigh,
                       cellScaleFactorNewMidLow, cellScaleFactorNewLow)
  proportionsVector <- c(predictedProportionsHigh, predictedProportionsMidHigh,
                         predictedProportionsMidLow, predictedProportionsLow)
  biasVector <- c(biasHigh, biasMidHigh, biasMidLow, biasLow)
  errorVector <- c(errorHigh, errorMidHigh, errorMidLow, errorLow)
  
  dfp2 <- do.call(rbind, lapply(seq(4), function(index){
    facetValue <- facetVector[index]
    newScalesValue <- newScalesVector[index]
    predictedProportion <- proportionsVector[index]
    biasValue <- biasVector[index]
    errorValue <- errorVector[index]
    
    cellScaleFactorChange <- newScalesValue-cellScaleFactorStart
    proportionChange <- predictedProportion-predictedProportionsStart
    biasChange <- biasValue-biasStart
    errorChange <- errorValue-errorStart
    markerExpressionChange <- 
      (markerExpression/newScalesValue)-markerExpression
    
    dfpNew <- data.frame(
      variable = 
        c("cellScaleFactor", "markerExpression", 
          "predictedProportion", "bias", "error"),
      value = 
        c(cellScaleFactorChange, markerExpressionChange, 
          proportionChange, biasChange, errorChange)
    )
    dfpNew$type <- facetValue
    return(dfpNew)
  })) |> as.data.frame()
  dfp2$Change <- ifelse(dfp2$value > 0, "Increase", "Decrease")
  dfp2$type <- factor(dfp2$type, levels = changeLevelsVector)
  dfp2$variable <- factor(dfp2$variable, 
                          levels = c("cellScaleFactor", "markerExpression",
                                     "predictedProportion", "bias", "error"))
  
  plot2 <- ggplot(dfp2, aes(x = variable, y = value, fill = Change)) + 
    geom_bar(stat="identity", color = "black") + theme_bw() +
    ylab(yLabelStringChanges) + facet_wrap(~type, nrow = 1) + 
    geom_hline(yintercept = 0) +
    #theme(axis.title.x = element_blank(),
    #      axis.text.x = element_text(angle=45,hjust=1),
    #      axis.ticks.y = element_blank(),
    #      axis.text.y = element_blank()) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    ggtitle(plotTitleStringChanges) +
    scale_fill_manual(breaks = changeLevelsVector, 
                      values=c("dodgerblue", "gold"))
  
  
  # return
  dataList <- list(
    plotDataValues = dfp,
    plotDataChanges = dfp2
  )
  plotList <- list(
    plotBarplotValues = plot1,
    plotBarplotChanges = plot2
  )
  returnList <- list(data = dataList,
                     plots = plotList)
  return(returnList)
}

listExample <- point_value_example()

listExample$plots$plotBarplotChanges
