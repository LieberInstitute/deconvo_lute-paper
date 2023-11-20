#!/usr/bin/env R

# Author: Sean Maden
#
# Performs deconvolution with the ABIS-seq signature matrix reference.
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
colnames(bulkExpression) <- gsub("_.*", "", colnames(bulkExpression))
colnames(bulkSummarizedExperiment) <- 
  gsub("_.*", "", colnames(bulkSummarizedExperiment))

# filter genes
# map symbols
bulkSummarizedExperimentNew <- yMapMarkers(bulkSummarizedExperiment)
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





#---------------------------------
# get common cell type id mappings
#---------------------------------

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
df.plot.tall.s13 <- as.data.frame(rbind(prop.scaled, prop.unscaled))

# get transpose
dfPlotTranspose <- t(df.plot.tall.s13) %>% as.data.frame()

#---------------------
# prep df.proportions
#---------------------
df.proportions <- experimentData[["p.true"]]

# append to df.plot.tall.s13
df.plot.tall.s13$true.proportions <- "NA"
df.plot.tall.s13$sample.id.format <- gsub("_.*", "", df.plot.tall.s13$sample.id)
for(type in colnames(df.proportions)){
  for(sample.id in unique(df.plot.tall.s13$sample.id.format)){
    filterTall <- 
    df.plot.tall.s13[
      df.plot.tall.s13[,"sample.id.format"]==sample.id, type]$true.proportions <- 
      df.proportions[sample.id,type]
  }
}




#------------------------------
# aggregate results proportions
#------------------------------
# get means on overlapping ids
df.plot.tall.mean <- 
  df.plot.tall.s13 %>% group_by(type) %>% 
  summarise(across(everything(), mean)) %>% t() %>% as.data.frame()
colnames(df.plot.tall.mean) <- df.plot.tall.mean[1,]
df.plot.tall.mean <- df.plot.tall.mean[seq(2,nrow(df.plot.tall.mean)),]
head(df.plot.tall.mean)

# get sds
df.plot.tall.sd <- df.plot.tall.s13 %>% 
  group_by(type) %>% summarise(across(everything(), sd)) %>%
  t() %>% as.data.frame()
colnames(df.plot.tall.sd) <- df.plot.tall.sd[1,]
df.plot.tall.sd <- df.plot.tall.sd[seq(2,nrow(df.plot.tall.sd)),]
head(df.plot.tall.sd)

# combine
colnames(df.plot.tall.sd) <- paste0("sd.", colnames(df.plot.tall.sd))
colnames(df.plot.tall.mean) <- paste0("mean.", colnames(df.plot.tall.mean))
df.tall <- as.data.frame(cbind(df.plot.tall.mean, df.plot.tall.sd))

# append true
df.tall$true.prop.mean.flow.cyto <- 
  df.tall$true.prop.sd.flow.cyto <- 
  df.tall$true.prop.mean.mrna.yield <- 
  df.tall$true.prop.sd.mrna.yield <-
  df.tall$s.cell.size <- "NA"
for(type in df.map$p.true.id){
  filter.tall <- gsub("\\.", " ", rownames(df.tall)) == type
  df.tall[filter.tall,"true.prop.mean.flow.cyto"] <- 
    df.proportions[df.proportions[,1]==type,"flow.cyto.mean"]
  df.tall[filter.tall,"true.prop.sd.flow.cyto"] <- 
    df.proportions[df.proportions[,1]==type,"flow.cyto.sd"]
  df.tall[filter.tall,"true.prop.mean.mrna.yield"] <- 
    df.proportions[df.proportions[,1]==type,"mrna.yield.mean"]
  df.tall[filter.tall,"true.prop.sd.mrna.yield"] <- 
    df.proportions[df.proportions[,1]==type,"mrna.yield.sd"]
  df.tall[filter.tall,"s.cell.size"] <- 
    s.vector[gsub("\\.", " ", names(s.vector))==type]
}
df.tall <- df.tall[!rownames(df.tall)=="sample.id",]

# format columns as numeric
for(c in seq(ncol(df.tall))){df.tall[,c] <- as.numeric(df.tall[,c])}
df.tall$cell.type <- rownames(df.tall)
# format cell size
df.tall$s.cell.size.log <- log(df.tall$s.cell.size)

# get df.wide
df.wide <- rbind(
  data.frame(mean = df.tall$mean.unscaled, sd = df.tall$sd.unscaled, 
             type = rep("unscaled", nrow(df.tall)), cell.type = rownames(df.tall),
             true.prop.mean.flow.cyto = df.tall$true.prop.mean.flow.cyto,
             true.prop.mean.mrna.yield = df.tall$true.prop.mean.mrna.yield,
             true.prop.sd.flow.cyto = df.tall$true.prop.sd.flow.cyto,
             true.prop.sd.mrna.yield = df.tall$true.prop.sd.mrna.yield),
  data.frame(mean = df.tall$mean.scaled, sd = df.tall$sd.scaled, 
             type = rep("scaled", nrow(df.tall)), cell.type = rownames(df.tall),
             true.prop.mean.flow.cyto = df.tall$true.prop.mean.flow.cyto,
             true.prop.mean.mrna.yield = df.tall$true.prop.mean.mrna.yield,
             true.prop.sd.flow.cyto = df.tall$true.prop.sd.flow.cyto,
             true.prop.sd.mrna.yield = df.tall$true.prop.sd.mrna.yield)
)
df.wide <- as.data.frame(df.wide)
df.wide <- df.wide[!df.wide$cell.type=="sample.id",]

# format columns as numeric
for(c in c(1,2,5,6,7,8)){df.wide[,c] <- as.numeric(df.wide[,c])}

# append bias and error
df.wide$true.prop.mean.flow.cyto.format <- 
  df.wide$true.prop.mean.flow.cyto/sum(na.omit(df.wide$true.prop.mean.flow.cyto))
df.wide$bias.flow.cyto <- df.wide$true.prop.mean.flow.cyto.format-
  df.wide$mean
df.wide$error.flow.cyto <- abs(df.wide$bias.flow.cyto)







#-----
# save
#-----
save.image("./env/02_abisseq/01_abisseq_script.RData")





