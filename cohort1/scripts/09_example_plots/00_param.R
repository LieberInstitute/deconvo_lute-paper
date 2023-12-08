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
  # cellScaleFactorsStart : starting cell scale factor value.
  # cellScaleFactorOffTypeValue : value of the difference for new simulation.
  # trueProportionValue : value of the true proportion across simulations.
  # markerExpressionStart : starting marker expression.
  # cellScaleFactorNew : new value of the scale factor.
  # bulkExpressionValue : value of bulk gene expressions
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
  # valuesList : List of starting simulation values.
  # roundValue : Second argument to base::round().
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
