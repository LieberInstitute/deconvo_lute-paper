#!/usr/bin/env/R

# Author: Sean Maden
#
# Run K5 experiments on ABIS-seq data
#
#
#
#
#
#
#


libv <- c('tidyr')
sapply(libv, library, character.only = TRUE)

#-----
# load
#-----

load("./env/02_abisseq/01_abisseq_script.RData")
load("./env/02_abisseq/02_proportions_s13_script.RData")
source("./source/lute_experiment.R")
source("./source/cell_mappings_helpers.R")

#------------------
# map k5 cell types
#------------------

vectorCellTypeMap <- c("T", "B", "Dendrocyte", "Plasma", "Monocyte", "Neutrophil", "NK")
vectorCellTypeStart <- colnames(zref)

mappingsTable <- cellLabelMappings(
  vectorCellTypeMap, vectorCellTypeStart, returnType = "list")[[1]]

zrefMapped <- cellLabelMappings(
  vectorCellTypeMap, vectorCellTypeStart, zref, 
  returnType = "df", summaryOperation = "mean")

# tests -- mappings for reference
#
#
#
# test label Monocyte
mappedLabel <- "Monocyte"

# test gene FXYD6, label Monocyte
geneName <- "FXYD6"
zrefMapped[geneName,mappedLabel]==rowMeans(
  zref[,which(colnames(zref) %in% mappingsTable[mappingsTable[,2]==mappedLabel,1])]
)[geneName]

# test gene NRG1, label Monocyte
geneName <- "NRG1"
zrefMapped[geneName,mappedLabel]==rowMeans(
  zref[,which(colnames(zref) %in% mappingsTable[mappingsTable[,2]==mappedLabel,1])]
)[geneName]

#
# test label Plasma
mappedLabel <- "Plasma"
# test gene FXYD6, label Plasma
geneName <- "FXYD6"
zrefMapped[geneName,mappedLabel]==
  zref[,which(colnames(zref) %in% 
                mappingsTable[mappingsTable[,2]==mappedLabel,1])][
                  which(rownames(zref)==geneName)]
# test gene NRG1, label Plasma
geneName <- "NRG1"
zrefMapped[geneName,mappedLabel]==
  zref[,which(colnames(zref) %in% 
                mappingsTable[mappingsTable[,2]==mappedLabel,1])][
                  which(rownames(zref)==geneName)]

# parse true proportions
#
#
#
# format from supplement
trueProportionsStart <- fc.proportions
rownames(trueProportionsStart) <- sampleNames <- fc.proportions[,1]
trueProportionsStart <- trueProportionsStart[,c(2:ncol(trueProportionsStart))]
trueProportionsStart <- trueProportionsStart[
  ,colnames(trueProportionsStart) %in% colnames(zref)]
# get fractions <1
trueProportionsStart <- t(apply(
  trueProportionsStart, 1, function(ri){ri/sum(ri)}))
trueProportionsStart <- as.data.frame(trueProportionsStart)
rownames(trueProportionsStart) <- sampleNames
# map labels
trueproportionsMapped <- cellLabelMappings(
  vectorCellTypeMap, colnames(trueProportionsStart), 
  trueProportionsStart, returnType = "df", summaryOperation = "sum")
trueproportionsMapped <- trueproportionsMapped[
  ,!colnames(trueproportionsMapped)=="NA"]
# get fractions <1
trueproportionsMapped <- t(apply(
  trueproportionsMapped, 1, function(ri){ri/sum(ri)}))
trueproportionsMapped <- as.data.frame(trueproportionsMapped)
rownames(trueproportionsMapped) <- sampleNames

# tests -- proportions
#
#
# test original proportions table
length(which(rowSums(trueProportionsStart)==1))==
  nrow(trueProportionsStart)
#
# test new proportions table
length(which(round(rowSums(trueproportionsMapped),1)==1))==
  nrow(trueproportionsMapped)


# parse cell type sizes
#
#
#
#
cellFactorsOriginal <- matrix(df.tall$s.cell.size,nrow=1)
colnames(cellFactorsOriginal) <- df.tall$cell.type
cellFactorsOriginal <- as.data.frame(cellFactorsOriginal)
dfCellScaleFactors <- cellLabelMappings(
  vectorCellTypeMap=vectorCellTypeMap, 
  vectorCellTypeStart=cellFactorsOriginal$cellTypesOriginal, 
  dfToMap=cellFactorsOriginal, returnType = "df", 
  summaryOperation = "mean"
)
cellScaleFactors <- as.numeric(dfCellScaleFactors)
names(cellScaleFactors) <- colnames(dfCellScaleFactors)
cellScaleFactors <- cellScaleFactors[!names(cellScaleFactors)=="NA"]

# test
#
#




#---------------
# run experiment
#---------------

bulkExpression <- as.matrix(assays(se)[["tpm"]])
bulkExpression <- bulkExpression[
  rownames(bulkExpression) %in% rownames(zref),]
colnames(bulkExpression) <- gsub("_.*", "", colnames(bulkExpression))
referenceExpression <- zrefMapped
referenceExpression <- referenceExpression[
  rownames(referenceExpression) %in% rownames(bulkExpression),]
trueCellTypeProportions <- trueproportionsMapped
colnames(trueCellTypeProportions) <- 
  paste0(colnames(trueCellTypeProportions), ".true")
trueCellTypeProportions$sample.id <- rownames(trueCellTypeProportions)

# get experiment results
#experimentList <- newExperimentList(
#  referenceExpression=referenceExpression,
#  trueCellTypeProportions=trueCellTypeProportions,
#  cellScaleFactors=cellScaleFactors,
#  bulkExpression=bulkExpression,
#  trueCellTypeProportionsSource="Flow cytometry",
#  typemarkerAlgorithmName=NULL
#)
#experimentResults <- evaluateExperiment(
#  experimentList, TRUE
#)

#--------------------------
# get the expression tables
#--------------------------
tpmReference <- zref

log2TpmReference <- 
  apply(
    tpmReference, 2, 
    function(cellName){log2(cellName+1)}) %>% as.data.frame()

scaleTpmReference <- scale(tpmReference) %>% as.data.frame()

scaleLog2TpmZref <- scale(log2TpmReference) %>% as.data.frame()

#-----
# save
#-----
save.image(file="./env/03_k5/01_run_script.RData")
