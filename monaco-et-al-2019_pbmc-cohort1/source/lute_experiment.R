#!/usr/bin/env R

# Author: Sean Maden
#
# Experiment functions. To be ported to object classes, eventually.
#
#

#' getRandomExperimentResults
#'
#' Get experiment results from using randomized experiment data.
#'
#'
getRandomExperimentResults <- function(
    sampleCount=10, typeCount=3, expressionMean=4, 
    markerPerType=50, cellsPerType=2000, cellScaleFactors=NULL){
  if(is(cellScaleFactors, "NULL")){
    cellScaleFactors <- rep(3, typeCount)
    cellScaleFactors[length(cellScaleFactors)] <- 10
    names(cellScaleFactors) <- paste0("type",seq(typeCount))
  } else{}
  
  # prep random data
  markerCount <- markerPerType*typeCount
  cellCount <- cellsPerType*typeCount
  singleCellExperiment <- randomSingleCellExperiment(
    numberGenes=markerCount, numberTypes=typeCount, numberCells=cellCount, 
    expressionMean=expressionMean)
  singleCellExperiment <- 
    singleCellExperiment[,sample(ncol(singleCellExperiment))]
  singleCellExperiment[["sample.id"]] <- 
    unlist(lapply(seq(sampleCount), function(i){
      rep(paste0("sample",i),cellCount/sampleCount)}))
  trueCellTypeProportions <- proportionsFromSingleCell(singleCellExperiment)
  trueCellTypeProportionsSource <- "SingleCellExperiment cell fractions"
  referenceExpression <- 
    referenceFromSingleCellExperiment(singleCellExperiment)
  uniqueSamples <- unique(singleCellExperiment[["sample.id"]])
  bulkExpression <- 
    do.call(cbind, lapply(uniqueSamples, function(sampleIter){
      filterSingleCell <- singleCellExperiment$sample.id==sampleIter
      ypb_from_sce(
        singleCellExperiment[,filterSingleCell], 
        assayName="counts", cellTypeVariable="celltype", 
        cellScaleFactors=cellScaleFactors
      )
    }))
  colnames(bulkExpression) <- uniqueSamples
  bulkExpression <- as.matrix(bulkExpression)
  experimentType <- "cellScaleFactorRescale"
  
  experimentListRandom <- newExperimentList(
    singleCellExperiment=NULL, referenceExpression=referenceExpression,
    bulkExpression=bulkExpression, cellScaleFactors=cellScaleFactors,
    experimentType=experimentType,
    trueCellTypeProportions=trueCellTypeProportions,
    trueCellTypeProportionsSource=trueCellTypeProportionsSource,
    assayName="counts", cellTypeVariable="celltype", 
    deconvolutionAlgorithmName="nnls", typemarkerAlgorithmName=NULL)
  return(
    evaluateExperiment(experimentListRandom)
  )
}

#' proportionsFromSingleCell
#'
#' Gets proportions, to be used as true cell type proportions, from 
#' SingleCellExperiment.
#'
#'
#'
proportionsFromSingleCell <- function(singleCellExperiment,
                                      sampleIdVariable="sample.id",
                                      cellTypeVariable="celltype"){
  
  sampleIdVector <- unique(singleCellExperiment[[sampleIdVariable]])
  uniqueCellTypes <- unique(singleCellExperiment[[cellTypeVariable]])
  
  # list available proportions
  proportionsList <- lapply(
    sampleIdVector, function(sample.id){
      filterSce <- 
        singleCellExperiment[[sampleIdVariable]]==sample.id
      sceFiltered <- singleCellExperiment[,filterSce]
      sampleProportions <- as.data.frame(
        prop.table(table(sceFiltered[[cellTypeVariable]])))
      sampleProportions$sample.id <- sample.id
      return(sampleProportions)
    })
  names(proportionsList) <- sampleIdVector
  
  # begin return table
  proportionsTable <- data.frame(sample.id = sampleIdVector)
  for(uniqueType in uniqueCellTypes){
    proportionsTable[,ncol(proportionsTable)+1] <- NA
    colnames(proportionsTable)[ncol(proportionsTable)] <- uniqueType
  }
  # append results
  for(sampleId in sampleIdVector){
    dataIter <- proportionsList[[sampleId]]
    for(uniqueType in uniqueCellTypes){
      proportionsTable[proportionsTable[,1]==sampleId,uniqueType] <- ifelse(
        uniqueType %in% dataIter[,1], dataIter[dataIter[,1]==uniqueType,2], 0
      )
    }
  }
  # format columns
  colnames(proportionsTable)[seq(2,ncol(proportionsTable))] <- 
    paste0(
      colnames(
        proportionsTable)[seq(2,ncol(proportionsTable))], ".true")
  
  return(
    proportionsTable
  )
}

#' newExperimentList
#'
#'
#'
#' @examples
#'
#'
newExperimentList <- function(singleCellExperiment=NULL,
                              referenceExpression=NULL,
                              bulkExpression=NULL,
                              cellScaleFactors=NULL,
                              trueCellTypeProportions=NULL,
                              trueCellTypeProportionsSource=NULL,
                              experimentType="cellScaleFactorRescale",
                              assayName="counts",
                              sampleIdVariable="sample.id",
                              cellTypeVariable="celltype",
                              deconvolutionAlgorithmName="nnls",
                              typemarkerAlgorithmName="findmarkers"){
  if(is(trueCellTypeProportions, "NULL")){
    trueCellTypeProportions <- proportionsFromSingleCell(
      singleCellExperiment, sampleIdVariable=sampleIdVariable,
      cellTypeVariable=cellTypeVariable)
  } else{}
  returnList <- list(
    singleCellExperiment=singleCellExperiment,
    referenceExpression=referenceExpression,
    bulkExpression=bulkExpression,
    cellScaleFactors=cellScaleFactors,
    trueCellTypeProportions=trueCellTypeProportions,
    trueCellTypeProportionsSource=trueCellTypeProportionsSource,
    experimentType=experimentType,
    assayName=assayName,
    cellTypeVariable=cellTypeVariable,
    deconvolutionAlgorithmName=deconvolutionAlgorithmName,
    typemarkerAlgorithmName=typemarkerAlgorithmName
  )
  return(returnList)
}

#' evaluateExperiment
#'
#'
#'
#'
evaluateExperiment <- function(experimentList, makePlots=TRUE){
  experimentType <- experimentList[["experimentType"]]
  if(experimentType=="cellScaleFactorRescale"){
    return(
      runExperimentCellScaleFactorsRescale(experimentList, makePlots=makePlots)
    )
  } else{
    stop("Error, experimentType not recognized.")
  }
  return(NULL)
}

#' runExperimentCellScaleFactorsRescale
#'
#'
#'
#'
runExperimentCellScaleFactorsRescale <- function(
    experimentList, makePlots = TRUE){
  # unpack
  cellScaleFactors <- experimentList[["cellScaleFactors"]]
  singleCellExperiment <- experimentList[["singleCellExperiment"]] 
  bulkExpression <- experimentList[["bulkExpression"]]
  cellScaleFactors <- experimentList[["cellScaleFactors"]]
  cellTypeVariable <- experimentList[["cellTypeVariable"]]
  assayName <- experimentList[["assayName"]]
  experimentType <- experimentList[["experimentType"]]
  trueProportions <- experimentList[["trueCellTypeProportions"]]
  deconvolutionAlgorithmName <- experimentList[["deconvolutionAlgorithmName"]]
  typemarkerAlgorithmName <- experimentList[["typemarkerAlgorithmName"]]
  
  # get conditions
  experimentConditionsList <- list(
    scaled = cellScaleFactors, unscaled = NULL
  )
  
  # get experiment results
  experimentResultsList <- lapply(
    seq(length(experimentConditionsList)), function(conditionIndex){
    conditionName <- names(experimentConditionsList)[conditionIndex]
    conditionValues <- experimentConditionsList[[conditionIndex]]
    if(is(singleCellExperiment, "NULL")){
      runResults <- lute(
        referenceExpression = referenceExpression, 
        bulkExpression = bulkExpression, 
        cellScaleFactors = conditionValues, # MAIN CONDITION VARIABLE
        typemarkerAlgorithm = typemarkerAlgorithmName,
        deconvolutionAlgorithm=deconvolutionAlgorithmName,
        cellTypeVariable = cellTypeVariable,
        assayName = assayName
      )
    } else{
      runResults <- lute(
        singleCellExperiment = singleCellExperiment, 
        bulkExpression = bulkExpression, 
        cellScaleFactors = conditionValues, # MAIN CONDITION VARIABLE
        typemarkerAlgorithm = NULL,
        deconvolutionAlgorithm=deconvolutionAlgorithmName,
        cellTypeVariable = cellTypeVariable,
        assayName = assayName
      )
    }
    return(
      list(runResults = runResults,
           experimentConditionName = conditionName,
           experimentConditionValues = conditionValues)
    )
  })
  names(experimentResultsList) <- paste0(
    "runResults;condition:",names(experimentConditionsList))
  
  # get proportions across conditions df
  proportionsTable <- getProportionsResultTable(
    experimentResultsList, experimentList)
  rmseResultTable <- rmseFromProportions(
    proportionsTable
  )
  corResultTable <- corFromProportions(
    proportionsTable, 
    cellTypesVector = names(cellScaleFactors),
    predictedString = "predicted",
    trueString = "true",
    corMethod = "pearson"
  )
  if(makePlots){
    plotResultList <- plotsFromCellScaleResults(
      proportionsTable, cellScaleFactors,
      stringTrue = "true", stringPred = "predicted"
    )
  } else{
    plotResultList <- NULL
  }
  
  # return results
  returnList <- list(
    experimentList=experimentList,
    experimentConditionsList=experimentConditionsList,
    experimentResultsList=experimentResultsList,
    proportionsTable=proportionsTable,
    rmseResultTable=rmseResultTable,
    corResultTable=corResultTable,
    plotResultList=plotResultList,
    experimentType=experimentType
  )
  return(returnList)
}

#' plotsFromCellScaleResults
#'
#' Lists new ggplot2 plot objects from results.
#'
#'
#'
#'
plotsFromCellScaleResults <- function(proportionsAcrossConditions,
                                      cellScaleFactors,
                                      stringTrue = "true",
                                      stringPred = "predicted",
                                      stringError = "error.true.pred"){
  # tall plot table
  dfPlot <- proportionsAcrossConditions
  columnsTrue <- c(
    "sample.id", "condition", 
    colnames(dfPlot)[
      grepl(
        paste0("\\.",stringTrue,"$"), colnames(dfPlot))])
  columnsPred <- c(
    "sample.id", "condition", 
    colnames(dfPlot)[
      grepl(
        paste0("\\.",stringPred,"$"), colnames(dfPlot))])
  
  dfPlotTallIter <- data.frame(sample.id = dfPlot$sample.id,
                           condition = dfPlot$condition)
  dfPlotTall <- do.call(rbind, 
                        lapply(names(cellScaleFactors), 
                               function(cell.type){
    dfPlotTallIterType <- dfPlotTallIter
    dfPlotTallIterType$celltype <- cell.type
    regexTypeTrue <- grepl(
      paste0(cell.type,"\\.", stringTrue, "$"),colnames(dfPlot))
    regexTypePred <- grepl(
      paste0(cell.type,"\\.", stringPred, "$"),colnames(dfPlot))
    regexTypeError <- grepl(
      paste0(cell.type,"\\.", stringError, "$"),colnames(dfPlot))
    dfPlotTallIterType$true <- dfPlot[,regexTypeTrue]
    dfPlotTallIterType$pred <- dfPlot[,regexTypePred]
    dfPlotTallIterType$error <- dfPlot[,regexTypeError]
    return(dfPlotTallIterType)
  }))
  dfPlotTall <- as.data.frame(dfPlotTall)
  # append cell sizes
  cellTypeVector <- unique(dfPlotTall$celltype)
  dfPlotTall$s.cell.size <- 1
  for(type in cellTypeVector){
    filterDf <- dfPlotTall$condition=="scaled" & dfPlotTall$celltype==type
    dfPlotTall[filterDf,]$s.cell.size <- cellScaleFactors[[type]]
  }
  
  # list plots
  # scatterplots
  listScatterTypes <- lapply(cellTypeVector, function(cell.type){
    scatterPlotList(dfPlotTall[dfPlotTall$celltype==cell.type,])
  })
  names(listScatterTypes) <- cellTypeVector
  listScatterAll <- scatterPlotList(dfPlotTall)
  # boxplots with jitters
  listBoxplotJitter <- boxplotJitterPlotList(dfPlotTall)
  
  # return
  returnList <- list(
    dfPlotTall=dfPlotTall,listScatterTypes=listScatterTypes,
    listScatterAll=listScatterAll,listBoxplotJitter=listBoxplotJitter
  )
  return(returnList)
}

#' plotWideScatter
#'
#' Get proportions scatterplots from tall plot table.
#'
scatterPlotList <- function(dfPlotTall){
  # scatterplots
  plot1 <- ggplot(dfPlotTall, aes(x = true, y = pred, color = sample.id)) + 
    geom_point(alpha = 0.5, size = 4) + facet_wrap(~condition*celltype) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_abline(slope = 1, intercept = 0) +
    xlim(0, 1) + ylim(0, 1) + 
    xlab("True") + ylab("Predicted")
  plot2 <- ggplot(dfPlotTall, aes(x = true, y = pred, color = sample.id)) + 
    geom_point(alpha = 0.5, size = 4) + facet_wrap(~condition*celltype) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_abline(slope = 1, intercept = 0) +
    xlab("True") + ylab("Predicted")
  plot3 <- ggplot(dfPlotTall, 
                  aes(x = true, y = pred, size = s.cell.size)) + 
    geom_point(alpha = 0.5) + facet_wrap(~condition*celltype) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_abline(slope = 1, intercept = 0) +
    xlim(0, 1) + ylim(0, 1) + 
    xlab("True") + ylab("Predicted")
  plot4 <- ggplot(dfPlotTall, 
                  aes(x = true, y = pred, size = s.cell.size)) + 
    geom_point(alpha = 0.5) + facet_wrap(~condition*celltype) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_abline(slope = 1, intercept = 0) +
    xlab("True") + ylab("Predicted")
  return(list(
    plot1=plot1, plot2=plot2, plot3=plot3, plot4=plot4
  ))
}

#' plotTallBoxplotJitter
#'
#' Get error boxplots with jitters from tall plot table.
#'
boxplotJitterPlotList <- function(dfPlotTall){
  # boxplots
  plot1 <- ggplot(dfPlotTall, 
                  aes(x = celltype, y = error)) +
    geom_jitter(alpha = 0.5, size = 4) + 
    geom_boxplot(color = "cyan", alpha = 0) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(~condition)
  plot2 <- ggplot(dfPlotTall, 
                  aes(x = celltype, y = error, color = sample.id)) +
    geom_jitter(alpha = 0.5, size = 4) + 
    geom_boxplot(color = "cyan", alpha = 0) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(~condition)
  plot3 <- ggplot(dfPlotTall, 
                  aes(x = celltype, y = error, 
                      color = sample.id, size = s.cell.size)) +
    geom_jitter(alpha = 0.5) + 
    geom_boxplot(color = "cyan", alpha = 0, size = 0.5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(~condition)
  plot4 <- ggplot(dfPlotTall, 
                  aes(x = sample.id, y = error, 
                      color = sample.id, size = s.cell.size)) +
    geom_jitter(alpha = 0.5) + 
    geom_boxplot(color = "cyan", alpha = 0, size = 0.5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(~condition*celltype)
  return(list(
    plot1=plot1, plot2=plot2, plot3=plot3, plot4=plot4
  ))
}


#' rmseFromProportions
#'
#' Gets RMSE summaries across samples and cell types, from proportions results 
#' table.
#'
#'
rmseFromProportions <- function(
    proportionsTable, errorString = "error.true.pred"){
  
  propColumns <- colnames(proportionsTable)
  regexErrorColumns <- grepl(paste0(".*\\.", errorString, "$"), propColumns)
  errorColumns <- propColumns[regexErrorColumns]
  
  # by cell type (across sample ids)
  conditions <- unique(proportionsTable$condition)
  rmseTypeTable <- do.call(rbind,
                          lapply(
                            conditions, function(condition){
    filteredTable <- proportionsTable[
      proportionsTable$condition==condition,]
    rmseType <- 
      apply(filteredTable[,errorColumns],2,rmseCall)
    rmseTypeTableIteration <- data.frame(
      sample.id.vector = paste0(
        unique(filteredTable$sample.id), 
        collapse = ";"),
      cell.type = gsub("\\..*", "", errorColumns),
      rmseType = rmseType
    )
    rmseTypeTableIteration$condition <- condition
    return(rmseTypeTableIteration)
  }))
  rownames(rmseTypeTable) <- NULL
  rmseTypeTable <- as.data.frame(rmseTypeTable)
  
  # by sample id (across cell types)
  rmseSampleTable <- do.call(rbind,
                           lapply(
                             conditions, function(condition){
               filteredTable <- proportionsTable[
                 proportionsTable$condition==condition,]
               rmseSample <- 
                 apply(filteredTable[,errorColumns],1,rmseCall)
               rmseSampleTableIteration <- data.frame(
                 celltype.vector = paste0(
                   unique(gsub("\\..*", "", errorColumns)),
                   collapse = ";"),
                 sample.id = filteredTable$sample.id,
                 rmseSample = rmseSample
               )
               rmseSampleTableIteration$condition <- condition
               return(rmseSampleTableIteration)
             }))
  rownames(rmseTypeTable) <- NULL
  rmseTypeTable <- as.data.frame(rmseTypeTable)
  
  # return
  returnList <- list(
    rmseTypeTable=rmseTypeTable,rmseSampleTable=rmseSampleTable)
  return(
    returnList
  )
}

#' rmseCall
#'
#' Main RMSE function
#'
rmseCall <- function(error){
  sqrt(mean(error^2))
}

#' getProportionsResultTable
#'
#' Gets the proportions results table of predictions from listed experiment 
#' results, experimentResultsList.
#' 
#' 
#' 
#'
getProportionsResultTable <- function(
    experimentResultsList, experimentList){
  trueProportions <- experimentList[["trueCellTypeProportions"]]
  cellScaleFactors <- experimentList[["cellScaleFactors"]]
  
  # bind proportions across conditions
  proportionsAcrossConditions <- do.call(rbind, lapply(
    seq(length(experimentResultsList)), function(index){
      resultsIteration <- experimentResultsList[[index]]
      conditionName <- resultsIteration$experimentConditionName
      dfProportionsIteration <- 
        resultsIteration$runResults[[
          "deconvolutionResults"
          ]]@predictionsTable
      dfProportions <- data.frame(
        sample.id = rownames(dfProportionsIteration),
        condition = rep(conditionName, nrow(dfProportionsIteration)))
      # parse types, allowing for missing types
      for(cellType in names(cellScaleFactors)){
        dfProportions[,ncol(dfProportions)+1] <- 0
        colnames(dfProportions)[ncol(dfProportions)] <- cellType
        if(cellType %in% colnames(dfProportionsIteration)){
          dfProportions[,ncol(dfProportions)] <- 
            dfProportionsIteration[,cellType]
        } else{}
      }
      return(dfProportions)
    }))
  proportionsAcrossConditions <- as.data.frame(proportionsAcrossConditions)
  filterColumns <- colnames(proportionsAcrossConditions) %in% 
    names(cellScaleFactors)
  colnames(proportionsAcrossConditions)[filterColumns] <-
    paste0(colnames(proportionsAcrossConditions)[filterColumns], ".predicted")
  # append true cell types proportions
  proportionsAcrossConditions <- 
    merge(
      proportionsAcrossConditions, trueProportions, value.name = "sample.id")
  
  # append bias and error
  proportionsAcrossConditions <- proportionsTableBiasAndError(
    proportionsAcrossConditions, names(cellScaleFactors)
  )
  return(
    proportionsAcrossConditions
  )
}

#' proportionsTableBiasAndError
#'
#' Adds bias and error columns to the proportions results table 
#' proportionsAcrossConditions.
#'
#'
#'
#'
proportionsTableBiasAndError <- function(proportionsAcrossConditions,
                                         cellTypesVector,
                                         predictedString = "predicted",
                                         trueString = "true"){
  proportionsTable <- proportionsAcrossConditions
  propColnames <- colnames(proportionsTable)
  regexTrueColumn <- grepl(paste0("\\.",trueString,"$"), propColnames)
  regexPredColumn <- grepl(paste0("\\.",predictedString,"$"), propColnames)
  for(cellType in cellTypesVector){
    regexTypeColumn <- grepl(paste0(cellType, "\\..*"), propColnames)
    regexTypeTrueColumn <- which(regexTrueColumn & regexTypeColumn)
    regexTypePredColumn <- which(regexPredColumn & regexTypeColumn)
    condition <- length(regexTypeTrueColumn)==1 & 
      length(regexTypePredColumn)==1
    if(condition){
      # bias
      proportionsTable[,ncol(proportionsTable)+1] <- "NA"
      colnames(
        proportionsTable)[
          ncol(proportionsTable)] <- 
        paste0(cellType,".bias.true.pred")
      proportionsTable[,paste0(cellType,".bias.true.pred")] <-
        proportionsTable[,regexTypeTrueColumn]-
        proportionsTable[,regexTypePredColumn]
      # error
      proportionsTable[,ncol(proportionsTable)+1] <- "NA"
      colnames(
        proportionsTable)[
          ncol(proportionsTable)] <- 
        paste0(cellType,".error.true.pred")
      proportionsTable[,paste0(cellType,".error.true.pred")] <-
        abs(proportionsTable[,paste0(cellType,".bias.true.pred")])
    }
  }
  return(proportionsTable)
}

#' corFromProportions
#'
#' Gets the correlation test results across experiment conditions. Calls 
#' corResultsCellType.
#'
#'
corFromProportions <- function(proportionsAcrossConditions, 
                               cellTypesVector,
                               predictedString = "predicted",
                               trueString = "true",
                               corMethod = "pearson"){
  if(length(unique(proportionsAcrossConditions$sample.id))<4){
    message("Warning, not enough samples. Skipping correlation tests.")
    return(NULL)
  } else{
    conditionsVector <- unique(proportionsAcrossConditions$condition)
    corResultsTable <- do.call(rbind, 
                               lapply(conditionsVector, function(condition){
      filterProportions <- proportionsAcrossConditions[
        proportionsAcrossConditions$condition==condition,]
      corResultsCellType(
        filterProportions, cellTypesVector, corMethod)
    }))
  }
  return(corResultsTable)
}

#' corResultsCellType
#'
#' Parses correlation test for an experiment condition. 
#'
corResultsCellType <- function(filterProportions, cellTypesVector,
                               predictedString = "predicted", 
                               trueString = "true", corMethod = "pearson"){
  propColnames <- colnames(filterProportions)
  regexTrueColumn <- grepl(paste0("\\.",trueString,"$"), propColnames)
  regexPredColumn <- grepl(paste0("\\.",predictedString,"$"), propColnames)
  resultTable <- do.call(rbind, 
            lapply(cellTypesVector, function(cellType){
              returnVector <- c(
                unique(filterProportions$condition), cellType, 
                paste0(filterProportions$sample.id, collapse = ';'),
                corMethod
              )
              regexTypeColumn <- grepl(paste0(cellType, "\\..*"), propColnames)
              regexTypeTrueColumn <- which(regexTrueColumn & regexTypeColumn)
              regexTypePredColumn <- which(regexPredColumn & regexTypeColumn)
              conditionPass <- length(regexTypeTrueColumn)==1 & 
                length(regexTypePredColumn)==1
              if(conditionPass){
                # corr test result
                corResult <- cor.test(filterProportions[,regexTypeTrueColumn], 
                                      filterProportions[,regexTypePredColumn], 
                                      method = corMethod)
                returnVector <- 
                  c(returnVector, corResult$estimate, corResult$p.value)
              } else{
                returnVector <- c(returnVector, "NA", "NA")
              }
              return(
                returnVector
              )
  }))
  colnames(resultTable) <- c(
    "condition", "cell.type", "sample.id.vector", 
    "method", "estimate", "pvalue")
  return(resultTable)
}

#'
#'
#'
#'
#'
postprocessExperiment <- function(resultsListPreprocess){
  
}