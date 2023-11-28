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

singleValueTestVariables <- function(cellScaleFactorsStart = 0.5,
                                     cellScaleFactorOffTypeValue = 1,
                                     trueProportionValue = 0.8,
                                     markerExpressionStart = 0.5,
                                     cellScaleFactorNew = 1,
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

parseExampleStartValues <- function(valuesList, roundValue = 2){
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
  matrixValues <- c(valuesList[["markerExpressionStart"]], 0.1, 0.1, 
                    valuesList[["markerExpressionStart"]])
  zrefExample <- matrix(matrixValues, nrow = 2)
  colnames(zrefExample) <- c("type1", "type2")
  zrefExampleScaled <- lute:::.zstransform(zrefExample, cellScaleFactorsStart)
  #bulkExpressionExample <- matrix(
  #  rep(valuesList[["bulkExpressionValue"]], 2), ncol = 1)
  bulkExpressionExample <- t(t(
    c(valuesList[["trueProportionValue"]], 
      1-valuesList[["trueProportionValue"]])) %*% 
      t(zrefExampleScaled))
  rownames(zrefExample) <- rownames(zrefExampleScaled) <- 
    rownames(bulkExpressionExample) <- paste0("gene", seq(nrow(zrefExample)))
  
  newParamStart <- 
    nnlsParam(bulkExpressionExample, zrefExample, 
              cellScaleFactorsStart) |>
    deconvolution()
  newParamNull <- 
    nnlsParam(bulkExpressionExample, zrefExample, 
              cellScaleFactorsNull) |>
    deconvolution()
  predictedProportionsStart <- newParamStart@predictionsTable[[1]]
  predictedProportionsNull <- newParamNull@predictionsTable[[1]]
  biasStart <- valuesList[["trueProportionValue"]]-predictedProportionsStart
  biasNull <- valuesList[["trueProportionValue"]]-predictedProportionsNull
  errorStart <- abs(biasStart)
  errorNull <- abs(biasNull)
  
  # make barplot
  dfp <- data.frame(
    cellScaleFactor = valuesList[["cellScaleFactorsStart"]],
    markerExpression = valuesList[["markerExpressionStart"]],
    trueProportion = valuesList[["trueProportionValue"]],
    predictedProportion = predictedProportionsStart,
    biasValue = biasStart,
    errorValue = errorStart
  )
  dfp <- melt(dfp)
  dfp$value <- round(dfp$value, roundValue)
  plot1 <- ggplot(dfp, aes(x = variable, y = value)) + 
    geom_bar(stat="identity", color = "black", fill = 'gray') + 
    theme_bw() + geom_hline(yintercept = 0) +
    theme(axis.text.x = element_text(angle=45,hjust=1),
          axis.title.x = element_blank()) +
    geom_text(aes(label = value), vjust = -1.5) +
    ylab("Value") + ylim(min(dfp$value-1), max(dfp$value+1)) +
    ggtitle("Example starting values")
  
  # return
  valuesList[["deconvoResultsStart"]] <- newParamStart
  valuesList[["bulkExpressionExample"]] <- bulkExpressionExample
  valuesList[["markerExpressionStartScaled"]] <- zrefExample[1,1]
  valuesList[["zrefExample"]] <- zrefExample
  valuesList[["zrefExampleScaled"]] <- zrefExampleScaled
  valuesList[["predictedProportionsStart"]] <- predictedProportionsStart
  valuesList[["predictedProportionsNull"]] <- predictedProportionsNull
  valuesList[["biasStart"]] <- biasStart
  valuesList[["biasNull"]] <- biasNull
  valuesList[["errorStart"]] <- errorStart
  valuesList[["errorNull"]] <- errorNull
  valuesList[["ggBarplotStart"]] <- plot1
  return(valuesList)
}

singleValueExample <- function(valuesList, conditionLabel = "", 
                               plotTitleString = ""){
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
      valuesList[["bulkExpressionExample"]], 
      valuesList[["zrefExample"]], 
      cellScaleFactorsNew) |>
    deconvolution()
  predictedProportionsNew <- newParamNew@predictionsTable[[1]]
  biasNew <- valuesList[["trueProportionValue"]]-predictedProportionsNew
  errorNew <- abs(biasNew)
  zrefNew <- 
    lute:::.zstransform(valuesList[["zrefExample"]], 
                        cellScaleFactorsNew)
  
  # get changes
  cellScaleFactorChange <- valuesList[["cellScaleFactorNew"]]-
    valuesList[["cellScaleFactorsStart"]]
  proportionChange <- predictedProportionsNew-
    valuesList[["predictedProportionsStart"]]
  biasChange <- biasNew-valuesList[["biasStart"]]
  errorChange <- errorNew-valuesList[["errorStart"]]
  markerExpressionChange <- zrefNew[1,1]-valuesList[["zrefExample"]][1,1]
  
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
    ggtitle(plotTitleString) +
    scale_fill_manual(breaks = c("Increase", "Decrease"), 
                      values=c("dodgerblue", "gold"))
  
  returnList <- list(
    valuesList = valuesList,
    deconvoResult = newParamNew,
    zrefNew = zrefNew,
    plotData = dfpNew,
    ggBarplotChange = plot2
  )
  
  return(returnList)
  
}

multiPanelPlots <- function(cellScaleFactorOffTypeValue = 10,
                            markerExpressionStart = 0.5,
                            cellScaleFactorsStart = 2,
                            trueProportionValue = 0.6,
                            cellScaleFactorVector = 
                              c(0.5, 1.5, 2.5, 3.5, 1),
                            labelVector = 
                              c("Decrease", "Slight Decrease", 
                                "Slight Increase", "Increase", "NULL")){
  # multiPanelPlots
  #
  # Get multiple plot panels from vector of changed cell types.
  #
  #
  #
  #
  
  listResults <- lapply(seq(length(cellScaleFactorVector)), function(index){
    cellScaleFactorIndex <- cellScaleFactorVector[index]
    labelIndex <- labelVector[index]
    valuesList <- singleValueTestVariables(
      cellScaleFactorsStart = cellScaleFactorsStart, 
      markerExpressionStart = markerExpressionStart,
      cellScaleFactorNew = cellScaleFactorIndex, 
      trueProportionValue = trueProportionValue)
    exampleResult <- singleValueExample(
      valuesList, paste0(labelIndex,"\ncellScaleFactor = ",cellScaleFactorIndex))
    
    return(
      list(result = exampleResult,
           plot = exampleResult$ggBarplotChange)
    )
  })
  names(listResults) <- labelVector
  
  # get plots formatted for grid arrange
  listPlots <- lapply(listResults, function(item){item[["plot"]]})
  names(listPlots) <- labelVector
  
  dfPlotAll <- do.call(rbind, lapply(listResults, function(item){
    item$result$plotData
  })) |> as.data.frame()
  dfPlotAll$conditionLabel <- 
    factor(dfPlotAll$conditionLabel, levels = unique(dfPlotAll$conditionLabel))
  ggMultiPanel <- ggplot(dfPlotAll, aes(x = variable, y = value, fill = Change)) + 
    geom_bar(stat="identity", color = "black") + theme_bw() +
    ylab("Change (New - Old)") + facet_wrap(~conditionLabel, nrow = 1) + 
    geom_hline(yintercept = 0) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    ggtitle("Affect of scale change") +
    scale_fill_manual(breaks = c("Increase", "Decrease"), 
                      values=c("dodgerblue", "gold"))
  
  # return
  returnList <- list(resultsList = listResults,
                     resultsPlotList = listPlots,
                     dfPlotAll = dfPlotAll,
                     ggMulti = ggMultiPanel)
  return(returnList)
}

#--------------------------------------------------

# multiple panels -- marker expression start is 0.5

#--------------------------------------------------

# get multiple panels
# view facet of multi panel plots
listMultiPlot <- multiPanelPlots()
listMultiPlot$ggMulti
barplotStart <- 
  listMultiPlot$resultsList$Decrease$result$valuesList$ggBarplotStart
grid.arrange(barplotStart, listMultiPlot$ggMulti, nrow = 1,
             layout_matrix = matrix(c(1,2,2,2),nrow=1))

# save
jpeg("./figures/09_example_plots/multipanel_example_barplots.jpg", width = 11, 
     height = 4, units = "in", res = 400)
grid.arrange(barplotStart, listMultiPlot$ggMulti, nrow = 1,
             layout_matrix = matrix(c(1,2,2,2),nrow=1))
dev.off()

# marker expression start is 1

# get multiple panels
# view facet of multi panel plots
listMultiPlot <- multiPanelPlots(markerExpressionStart = 1)
barplotStart <- 
  listMultiPlot$resultsList$Decrease$result$valuesList$ggBarplotStart
grid.arrange(barplotStart, listMultiPlot$ggMulti, nrow = 1,
             layout_matrix = matrix(c(1,2,2,2),nrow=1))

#--------------

# single panels

#--------------
cellScaleFactorsStart <- 1
cellScaleFactorNew <- 3
trueProportionValue <- 0.2
valuesList <- singleValueTestVariables(
  cellScaleFactorsStart = cellScaleFactorsStart, 
  cellScaleFactorNew = cellScaleFactorNew, 
  trueProportionValue = trueProportionValue)
exampleResult <- singleValueExample(
  valuesList, paste0("Increase\ncellScaleFactor = ",cellScaleFactorNew))
grid.arrange(exampleResult$valuesList$ggBarplotStart,
             exampleResult$plot, nrow = 1)

cellScaleFactorsStart <- 1
cellScaleFactorNew <- 0.5
trueProportionValue <- 0.2
valuesList <- singleValueTestVariables(
  cellScaleFactorsStart = cellScaleFactorsStart, 
  cellScaleFactorNew = cellScaleFactorNew, 
  trueProportionValue = trueProportionValue)
exampleResult <- singleValueExample(
  valuesList, paste0("Decrease\ncellScaleFactor = ",cellScaleFactorNew))
grid.arrange(exampleResult$valuesList$ggBarplotStart,
             exampleResult$plot, nrow = 1)


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
