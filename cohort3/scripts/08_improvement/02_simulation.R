#!/usr/bin/env R

# Author: Sean Maden
#
# Pseudobulk simulation using data from ABIS-seq and flow cytometry (FC) true 
# proportions, and cell sizes from library sizes (column sums) of ABIS-seq reference.
#
#
#
#
#


libv <- c("lute", "dplyr")
sapply(libv, library, character.only = TRUE)

# load 
# contains experiment data, cell size scale factors
load("./env/06_top10markers/01_run_script.RData")

# run simulation

ptrue <- experimentData$p.true
cellLabelLarge <- "Plasmablasts"
cellLabelNotLarge <- "Non-plasmablast"
totalCells <- 10000
numberGenes <- nrow(refTpmFilter)
expressionMean <- mean(refTpmFilter)
cellScaleFactors <- c(cellSizes[cellLabelLarge],
                      median(cellSizes[!names(cellSizes)==cellLabelLarge]))
names(cellScaleFactors) <- c(cellLabelLarge, cellLabelNotLarge)

listPseudoBulk <- lapply(seq(ncol(ptrue)), function(index){
  isLargeCell <- rownames(ptrue)==cellLabelLarge
  # cell labels from fractions
  fractLargeCell <- ptrue[isLargeCell,index]
  totalLargeCells <- round(fractLargeCell*totalCells, 0)
  totalNotLargeCells <- round((1-fractLargeCell)*totalCells, 0)
  newCellLabels <- c(rep(cellLabelLarge, totalLargeCells),
                     rep(cellLabelNotLarge, totalNotLargeCells))
  # fraction vector
  #pFractionVector <- c(fractLargeCell, 1-fractLargeCell)
  #names(pFractionVector) <- c(cellLabelLarge, cellLabelNotLarge)
  pFractionVector <- c(totalLargeCells/totalCells, totalNotLargeCells/totalCells)
  names(pFractionVector) <- c(cellLabelLarge, cellLabelNotLarge)
  # randomize sce data
  scePseudo <- randomSingleCellExperiment(
    numberGenes = numberGenes,
    numberCells = totalCells,
    numberType = 2,
    fractionTypes = pFractionVector,
    expressionMean = expressionMean
  )
  scePseudo[,scePseudo$celltype=="type1"]$celltype <- cellLabelLarge
  scePseudo[,scePseudo$celltype=="type2"]$celltype <- cellLabelNotLarge
  yPseudo <- lute::ypb_from_sce(
    scePseudo, cellScaleFactors = cellScaleFactors)
  
})

fractPlasmablast <- 
numberCells <- 1000

expressionMean <- mean(refTpm)
sc <- randomSingleCellExperiment(
  numberType = 2,
  expressionMean = expressionMean,
  numberCells = 1000
)





