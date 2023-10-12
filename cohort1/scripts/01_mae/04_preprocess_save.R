#!/usr/bin/env R

# Author: Sean Maden
# 
# Preprocess MultiAssayExperiment by applying filters and saving MAE set for analysis.
#
#
#
#
#

filter.rnascope.confidence <- "Low"

#-----
# load
#-----
# mae input
mae.in.path <- "./outputs/01_mae/mae_allsamples_append.rda"
mae <- get(load(mae.in.path))
# test mae input
num.samples.mae.input <- 22
length(unique(colData(mae)$sample.id))==num.samples.mae.input

# rnascope confidence annotations
cd.id <- get(
  load(
    "./outputs/01_mae/sample_qc_df.rda"))
samples.to.remove <- unique(cd.id[cd.id$remove.low==TRUE,]$sample.id)

# samples to keep
sample.vector.keep <- colData(mae)$sample.id
sample.vector.keep <- sample.vector.keep[!sample.vector.keep %in% samples.to.remove]

# test filtered samples
num.samples.keep <- 19
num.samples.filter <- 3
length(samples.to.remove)==num.samples.filter
length(sample.vector.keep)==num.samples.keep







#---------------------------------
# 1. defines the sample ids filter
#---------------------------------
total.source.ids.to.keep <- 13

# sample id vector to keep
sample.vector.keep <- unique(
  na.omit(cd.id[!cd.id$sample.id %in% samples.to.remove,]$sample.id))

sample.vector.keep <- intersect(
  unique(mae[[1]]$Sample),
  unique(sample.vector.keep)
) # unique with k2 sce samples

#sample.vector.keep <- unique(
#  unique(mae[[1]]$Sample),
#  sample.vector.keep
#) # unique with k2 sce samples

length(sample.vector.keep)==total.source.ids.to.keep














#--------------------------------------
# 2. remove all non-overlapping samples
#--------------------------------------
mae.old <- mae

sample.mae.map.filter.final <- 
  colData(mae)$sample.id %in% sample.vector.keep

mae <- mae[,sample.mae.map.filter.final,]

# tests filtered mae
final.samples.count.mae <- 13
length(unique(colData(mae)$sample.id))==final.samples.count.mae

# tests filtered mae
names.test.mae1 <- c("snrnaseq.k2.all")
for(name1 in names.test.mae1){
  message(name1)
  message(length(unique(mae[[name1]]$Sample))==final.samples.count.mae)
}

names.test.mae2 <- c("snrnaseq.k3.all", 
                    "snrnaseq.k4.all",
                    "bulk.rnaseq", 
                    "bulk.rpkm.rnaseq",
                    "cell.sizes", 
                    "sce.img",
                    "bulk.pb.k2", "bulk.pb.k3", "bulk.pb.k4")
for(name2 in names.test.mae2){
  message(name2)
  message(
    length(
      unique(mae[[name2]]$Sample))==final.samples.count.mae)
}

names.test.mae3 <- c("cell.sizes")
for(name in names.test.mae3){
  message(name)
  message(length(unique(mae[[name]]["sample_id",]))==mae.final.samples.count)
}














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
