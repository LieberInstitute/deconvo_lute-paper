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


libv <- c("tidyr", "lute")
sapply(libv, library, character.only = TRUE)
source("./scripts/03_k5/00_param.R")

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

vectorCellTypeMap <- 
  c("T", "B", "Dendrocyte", "Plasma", "Monocyte", "Neutrophil", "NK")
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
#
#--------------------------

log2TpmReference <- 
  apply(
    tpmReference, 2, 
    function(cellName){log2(cellName+1)}) %>% as.data.frame()




#--------------------------
# get the expression tables
#--------------------------
tpmReference <- zref

# list quantile results
#
#

getQuantileTablesFromReferenceExpression(
  tpmReference, "TPM", 10, seq(0,1,0.1)
)[[3]]

getQuantileTablesFromReferenceExpression(
  scaleTpmReference, "Z TPM", 10, seq(0,1,0.1)
)[[3]]

getQuantileTablesFromReferenceExpression(
  log2TpmReference, "log2 TPM", 10, seq(0,1,0.1)
)[[3]]

getQuantileTablesFromReferenceExpression(
  scaleLog2TpmZref, "Z log2 TPM", 10, seq(0,1,0.1)
)[[3]]

scaleTpmReference <- scale(tpmReference) %>% as.data.frame()

scaleLog2TpmReference <- scale(log2TpmReference) %>% as.data.frame()

#
#
#
listQuantileTables <- 
  getQuantileTablesFromReferenceExpression(tpmReference)
# test
identical(
  ifelse(
    listQuantileTables$booleanTable[1,seq(10)], 1, 0),
  listQuantileTables$numericTable[1,seq(10)]
) # TRUE


#-----------------------
# get cluster results
#-----------------------
clusterTpmReference <- prcomp(tpmReference)
clusterTpmLog2Reference <- prcomp(log2TpmReference)
clusterScaleTpmReference <- prcomp(scaleTpmReference)
clusterScaleTpmLog2Reference <- prcomp(scaleLog2TpmZref)


#-------------------------
# get quantiles categories
#-------------------------

# quantiles for analysis scale (log2 [tpm + 1])
#
listQuantilesLog2TpmRef <- 
getQuantileTablesFromReferenceExpression(
  log2TpmReference, "log2 TPM", 10, seq(0,1,0.1)
)
# gets tall version for plots and
#
dfTallLog2TpmRef <- melt(
  listQuantilesLog2TpmRef[["booleanTable"]])


#--------------------------
# prepare data for heatmaps
#--------------------------
# map plot cols to cell types
# use 2 cell type mappings
#




mappingsTable$colorLabel1 <- mappingsTable$celltype1 %>%
  as.factor() %>% as.numeric()
mappingsTable$colorLabel2 <- mappingsTable$celltype2 %>%
  as.factor() %>% as.numeric()

# get markers vector with colors
markerColors <- data.frame(
  marker=rownames(scaleLog2TpmReference)
)
markerColors$cellType <- 
  markerColors$markerColor1 <- 
  markerColors$markerColor2 <- "NA"
for(type in colnames(scaleLog2TpmReference)){
  typeMarkerNames <- 
    dfTallLog2TpmRef[
      dfTallLog2TpmRef$Var2==type & dfTallLog2TpmRef$value==TRUE,]$Var1
  markerColors[markerColors$marker %in% typeMarkerNames,]$cellType1 <- type
  markerColors[markerColors$marker %in% typeMarkerNames,]$cellType1 <- 
    mappingsTable[mappingsTable$celltype1==type,]$celltype2
  markerColors[markerColors$marker %in% typeMarkerNames,]$markerColor1 <- 
    mappingsTable[mappingsTable$celltype1==type,]$colorLabel1
  markerColors[markerColors$marker %in% typeMarkerNames,]$markerColor2 <- 
    mappingsTable[mappingsTable$celltype1==type,]$colorLabel2
}

dfPlotHeatmapLog2TpmRef <- dfTallLog2TpmRef %>% 
  group_by(Var2) %>% 
  count(value) %>% 
  mutate(prop = prop.table(n))
colnames(dfPlotHeatmapLog2TpmRef) <- 
  c("cellType", "isTypeMarker", "markerCount", "markerProportion")










#-----
# save
#-----
save.image(file="./env/03_k5/01_run_script.RData")
