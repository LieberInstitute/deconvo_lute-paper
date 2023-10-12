#!/usr/bin/env R

#
# Preprocess MultiAssayExperiment
#
#
#
#
#


filter.rnascope.confidence <- "Low"

#-----
# load
#-----
mae.in.path <- "./outputs/01_mae/mae_allsamples_append.rda"
mae <- mae.all <- get(load(mae.in.path))
dim(colData(mae))

# rnascope confidence annotations
cd.id <- get(load("./outputs/01_mae/sample_qc_df.rda"))
samples.to.remove <- 
  cd.id$confidence.star == "Low" |
  cd.id$confidence.circle == "Low"
sample.vector.keep <- cd.id[!samples.to.remove,]$sample.id






#---------------------------------------------
# 1. filter on rnascope confidence annotations
#---------------------------------------------

# filter sce.img
sce.img <- mae[["sce.img"]]
dim(sce.img)
filter.sce <- colData(sce.img)$Sample %in% sample.vector.keep
sce.img <- sce.img[,filter.sce]
dim(sce.img)
mae[["sce.img"]] <- sce.img

# filter cell.sizes
cell.sizes <- mae[["cell.sizes"]]
dim(cell.sizes)
filter.cell.sizes <- cell.sizes["sample_id",] %in% sample.vector.keep
cell.sizes <- cell.sizes[,filter.cell.sizes]
dim(cell.sizes)
mae[["cell.sizes"]] <- cell.sizes


#--------------------------------------
# 2. remove all non-overlapping samples
#--------------------------------------
mae.old <- mae

sample.map.mae <- mae@sampleMap

overlapping.samples.3.platforms <- intersect(
  intersect(
    sample.map.mae[sample.map.mae$assay=="snrnaseq.k2.all",]$primary,
    sample.map.mae[sample.map.mae$assay=="cell.sizes",]$primary
  ),
  sample.map.mae[sample.map.mae$assay=="bulk.rnaseq",]$primary
)

#filter.samples.3platforms <- 
#  sample.map.mae$primary %in% overlapping.samples.3.platforms
#sample.map.mae.final <- sample.map.mae[filter.samples.3platforms,]

sample.mae.map.filter.final <- 
  colData(mae)$sample.id == unique(overlapping.samples.3.platforms)

mae <- mae[,sample.mae.map.filter.final,]

dim(colData(mae))










#----------------------------------------
# tests, 
# sample source id counts by mae platform
#----------------------------------------

# tests for rnascope input data
# test
number.samples.to.keep <- 13
length(sample.vector.keep)==number.samples.to.keep
# test: samples quantity post filter
expected.rnascope.samples.post.filter <- 16
sample.vector.keep <- cd.id[cd.id$remove.low==FALSE,]$sample.id
length(sample.vector.keep)==expected.rnascope.samples.post.filter



# tests
# test: number of samples after rnascope filter
length(unique(mae[["cell.sizes"]]["sample_id",]))==expected.rnascope.samples.post.filter
# test: number of snrnaseq and rnascope samples overlapping, post filter
expected.rnascope.snrnaseq.samples.post.filter <- 13
length(intersect(unique(mae[["snrnaseq.k2.all"]]$Sample),
                 unique(mae[["cell.sizes"]]["sample_id",])))==
  expected.rnascope.snrnaseq.samples.post.filter
length(intersect(
  sample.map.mae[sample.map.mae$assay=="snrnaseq.k2.all",]$primary,
  sample.map.mae[sample.map.mae$assay=="cell.sizes",]$primary
))==expected.rnascope.snrnaseq.samples.post.filter

# test: number of samples intersecting (BULK, RNASCOPE.FILTER, SNRNASEQ)
expected.rnascope.snrnaseq.samples.post.filter <- 13
length(
  intersect(
    intersect(
      unique(mae[["snrnaseq.k2.all"]]$Sample),
      unique(mae[["cell.sizes"]]["sample_id",])
    ),
    unique(mae[["bulk.rnaseq"]]$batch.id2)
  )
) == expected.rnascope.snrnaseq.samples.post.filter
length(
  intersect(
    intersect(
      sample.map.mae[sample.map.mae$assay=="snrnaseq.k2.all",]$primary,
      sample.map.mae[sample.map.mae$assay=="cell.sizes",]$primary
    ),
    sample.map.mae[sample.map.mae$assay=="bulk.rnaseq",]$primary
  )
) == expected.rnascope.snrnaseq.samples.post.filter







#-----
# save
#-----

dim(colData(mae))
mae.out.path <- "./outputs/01_mae/mae_analysis_append.rda"
save(mae, file = mae.out.path)
