#!/usr/bin/env R

# Author: Sean Maden
#
# Performs deconvolution with the ABIS-seq signature matrix reference.
#
#
#

libv <- c("lute")
sapply(libv, library, character.only = TRUE)

#-----
# load
#-----
zref <- read.csv("./data/zref_abisseq.csv")
ptrue <- read.csv("./data/ptrue_abis.csv")
load("./env/01_tpm_summaries/01_read_script.RData")

#--------------
# format inputs
#--------------
# zref
zref <- zref[!zref[,1]=="12:00 AM",]
rownames(zref) <- zref[,1]
zref <- zref[,c(2:ncol(zref))]

# se
rownames(se) <- rowData(se)$gene_symbol
filter.rows <- !(duplicated(rownames(se))|is.na(rownames(se)))
se <- se[filter.rows,]

# ptrue
df.proportions <- ptrue
df.proportions <- df.proportions[c(3:nrow(df.proportions)),]
colnames(df.proportions) <- c("cell.type", 
                              "flow.cyto.mean", 
                              "flow.cyto.sd",
                              "mrna.yield.mean",
                              "mrna.yield.sd")

# get id mappings
map.vector1 <- gsub("\\.", " ", colnames(prop.scaled))
map.vector2 <- df.proportions$cell.type
common.id <- intersect(map.vector1, map.vector2)
df.map <- data.frame(p.true.id = common.id,
                     map.vector2 = common.id)

#-----------------------
# get cell scale factors
#-----------------------
# marker library size
s.vector <- colSums(zref)

#------------
# run deconvo
#------------

result.unscaled <- lute(
  z = as.matrix(zref), 
  y = as.matrix(assays(se)[["tpm"]]), 
  assay.name = 'tpm',
  typemarker.algorithm = NULL
)

result.scaled <- lute(
  z = as.matrix(zref), 
  y = as.matrix(assays(se)[["tpm"]]),
  s = s.vector,
  assay.name = 'tpm',
  typemarker.algorithm = NULL
)

# make df.plot.tall
prop.unscaled <- result.unscaled[[1]]@predictions.table
prop.scaled <- result.scaled[[1]]@predictions.table

prop.unscaled$sample.id <- rownames(prop.unscaled)
prop.unscaled$type <- "unscaled"
prop.scaled$sample.id <- rownames(prop.scaled)
prop.scaled$type <- "scaled"




#-----
# save
#-----
save.image("./env/02_abisseq/01_abisseq_script.RData")
