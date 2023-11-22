#!/usr/bin/env R

# Author: Sean Maden
#
# Performs deconvolution with the ABIS-seq signature matrix reference.
#
#
#
#
#


libv <- c("lute", "dplyr")
sapply(libv, library, character.only = TRUE)

# load 
# 01_run_script
# contains experiment data, cell size scale factors
load("./env/06_top10markers/01_run_script.RData")
source("./source/00_read_experiment_data.R")



#--------------
# format inputs
#--------------
# cell type reference
referenceExpression <- as.matrix(experimentData[["z"]])

# get y/bulk expression -- only PBMC
bulkExpression <- as.matrix(experimentData[["y.table"]])
bulkExpression <- bulkExpression[,grepl("PBMC", colnames(bulkExpression))]
dim(bulkExpression)
bulkSummarizedExperiment <- experimentData[["y.se"]]
bulkSummarizedExperiment <- 
  bulkSummarizedExperiment[,grepl("PBMC", colnames(bulkSummarizedExperiment))]
dim(bulkSummarizedExperiment)

# format sample ids
# check colnames on table and se data
identical(colnames(bulkSummarizedExperiment), colnames(bulkExpression))
# replace underscore characters
colnames(bulkExpression) <- gsub("_.*", "", colnames(bulkExpression))
colnames(bulkSummarizedExperiment) <- 
  gsub("_.*", "", colnames(bulkSummarizedExperiment))
# replace x characters
colnames(bulkExpression) <- gsub("^X", "", colnames(bulkExpression))
colnames(bulkExpression)
# compare fc proportions ids
fcProportionsIds <- colnames(experimentData$p.true)
length(intersect(colnames(bulkExpression),fcProportionsIds))
# set new sampleids as colnames
colnames(bulkSummarizedExperiment) <- colnames(bulkExpression)
# subset on overlapping ids
sampleIdOverlap <- intersect(colnames(bulkSummarizedExperiment), fcProportionsIds)
bulkSummarizedExperiment <- bulkSummarizedExperiment[,sampleIdOverlap]
bulkExpression <- bulkExpression[,sampleIdOverlap]

# filter genes
# map symbols
bulkSummarizedExperimentNew <- yMapMarkers(bulkSummarizedExperiment)
save(bulkSummarizedExperimentNew, file = "./outputs/08_improvements/bulkWithGeneSymbols.rda")

# filter duplicates
filterDuplicatedGenes <- duplicated(
  rowData(bulkSummarizedExperimentNew)$gene_symbol)
bulkSummarizedExperimentNew <- 
  bulkSummarizedExperimentNew[!filterDuplicatedGenes]
# filter markers overlaps
rownames(bulkSummarizedExperimentNew) <- 
  rowData(bulkSummarizedExperimentNew)$gene_symbol
markersOverlaps <- 
  rownames(bulkSummarizedExperimentNew) %in% rownames(referenceExpression)
bulkSummarizedExperimentNew <- 
  bulkSummarizedExperimentNew[markersOverlaps,]

# inspect
dim(bulkSummarizedExperimentNew)
length(intersect(
  rownames(bulkSummarizedExperimentNew), rownames(referenceExpression)))
bulkExpressionNew <- as.matrix(assays(bulkSummarizedExperimentNew)[["tpm"]])





#------------
# run deconvo
#------------
result.unscaled <- lute(
  referenceExpression = referenceExpression, 
  bulkExpression = bulkExpressionNew,
  assayName = 'tpm',
  typemarkerAlgorithm = NULL
)



result.scaled <- lute(
  referenceExpression = referenceExpression, 
  bulkExpression = bulkExpressionNew,
  cellScaleFactors = cellSizes,
  assayName = 'tpm',
  typemarkerAlgorithm = NULL
)





#------------------
# make df.plot.tall
#------------------

prop.unscaled <- result.unscaled[["deconvolutionResults"]]@predictionsTable
prop.scaled <- result.scaled[["deconvolutionResults"]]@predictionsTable

# get common cell type id mappings
map.vector1 <- gsub("\\.", " ", colnames(prop.scaled))
map.vector2 <- gsub("\\.", " ", colnames(prop.unscaled)) # df.proportions$cell.type
common.id <- intersect(map.vector1, map.vector2)
df.map <- data.frame(
  p.true.id = common.id, map.vector2 = common.id)

# filter common ids
filter.columns <- gsub("\\.", " ", colnames(prop.scaled)) %in% df.map[,1]
prop.scaled <- prop.scaled[,filter.columns]
prop.unscaled <- prop.unscaled[,filter.columns]

# bind results
prop.unscaled$sample.id <- rownames(prop.unscaled)
prop.unscaled$type <- "unscaled"
prop.scaled$sample.id <- rownames(prop.scaled)
prop.scaled$type <- "scaled"

# get plot data
# get tall table (SAMPLES X CELL_TYPE)
dfPlotTall <- rbind(prop.scaled, prop.unscaled) %>% as.data.frame()
# get wide table (CELL_TYPE X SAMPLES)
dfPlotTranspose <- t(dfPlotTall) %>% as.data.frame()

# append proportions
dfProportions <- experimentData[["p.true"]]
# check samples overlap
length(intersect(colnames(dfProportions), dfPlotTall$sample.id))
samplesKeep <- intersect(colnames(dfProportions), unique(dfPlotTall$sample.id))

# get very tall plot of sample X proportion
dfPlotTallTall <- reshape2::melt(dfPlotTall)
colnames(dfPlotTallTall) <- 
  c("sample.id", "condition", "cellType", "proportion")
trueProportionVarname <- "trueProportion"
dfPlotTallTall[,ncol(dfPlotTallTall)+1] <- "NA"
colnames(dfPlotTallTall)[ncol(dfPlotTallTall)] <- trueProportionVarname
for(cellType in rownames(dfProportions)){
  for(sample.id in unique(dfPlotTallTall$sample.id)){
    filterTallTall <- 
      dfPlotTallTall[,"sample.id"]==sample.id & 
      dfPlotTallTall[,"cellType"]==cellType
    dfPlotTallTall[filterTallTall, trueProportionVarname] <- 
      dfProportions[cellType, sample.id]
  }
}
# filter samples
filterSamples <- dfPlotTallTall$sample.id %in% samplesKeep
dfPlotTallTall <- dfPlotTallTall[filterSamples,]
dfPlotTallTall$trueProportion <- as.numeric(dfPlotTallTall$trueProportion)

#-----------------
# save environment
#-----------------
save.image("./env/08_improvement/01_run_script.RData")

