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

point_value_example <- function(cellScaleFactorStart = 10,
                                cellScaleFactorOffTypeValue = 3,
                                trueProportionValue = 0.8,
                                markerExpression = 5,
                                cellScaleFactorNewHigh = 0.95,
                                cellScaleFactorNewMidHigh = 0.75,
                                cellScaleFactorNewMidLow = 0.25,
                                cellScaleFactorNewLow = 0.05,
                                changeNull = 0){
  # point_value_example
  #
  # Get example from simulation using provided point values
  #
  #
  #
  
  changeHigh <- cellScaleFactorNewHigh-cellScaleFactorStart
  changeMidHigh <- cellScaleFactorNewMidHigh-cellScaleFactorStart
  changeMidLow <- cellScaleFactorNewMidLow-cellScaleFactorStart
  changeLow <- cellScaleFactorNewLow-cellScaleFactorStart
  labelHigh <- paste0("Increase\n(cellScaleFactor = ",cellScaleFactorNewHigh,")")
  labelMidHigh <- paste0("Moderate Inc.\n(cellScaleFactor = ",cellScaleFactorNewMidHigh,")")
  labelMidLow <- paste0("Moderate Dec.\n(cellScaleFactor = ",cellScaleFactorNewMidLow,")")
  labelLow <- paste0("Decrease\n(cellScaleFactor = ",cellScaleFactorNewLow,")")
  
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
  cellScaleFactorsMidHigh <- c(cellScaleFactorNewMidHigh, cellScaleFactorOffTypeValue)
  cellScaleFactorsMidLow <- c(cellScaleFactorNewMidLow, cellScaleFactorOffTypeValue)
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
    predictedProportion = predictedProportionValue,
    biasValue = biasValue,
    errorValue = errorValue
  )
  dfp <- melt(dfp)
  
  plot1 <- ggplot(dfp, aes(x = variable, y = value)) + 
    geom_bar(stat="identity", color = "black", fill = 'gray') + 
    theme_bw() + ylim(-0.2, 0.75) + geom_hline(yintercept = 0) +
    theme(axis.text.x = element_text(angle=45,hjust=1),
          axis.title.x = element_blank()) +
    geom_text(aes(label = value), vjust = -1.5) +
    ylab("Value")
  
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
      (markerExpression*newScalesValue)-markerExpression
    
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
  
  plot2 <- ggplot(dfp2, aes(x = variable, y = value, fill = Change)) + 
    geom_bar(stat="identity", color = "black") + theme_bw() +
    ylab(yLabelString) + facet_wrap(~type, nrow = 1) + 
    geom_hline(yintercept = 0) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle=45,hjust=1),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank()) +
    ggtitle(plotTitleString) +
    scale_fill_manual(breaks = changeLevelsVector, 
                      values=c("dodgerblue", "gold")) +
    ylim(-1,1)
  
  
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
                     plots = plotsList)
  return(returnList)
}


listExample <- point_value_example()



