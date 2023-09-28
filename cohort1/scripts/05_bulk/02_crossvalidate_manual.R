#!/usr/bin/env R

# Author: Sean Maden
#
# Cross-validate normalization tests. Use best normalization procedure on bulk validation samples subset
#
#

libv <- c("scuttle")
sapply(libv, library, character.only = TRUE)

#--------------------
# get validation data
#--------------------
# Prepares MultiAssayExperiment object
sample.id.validate <- get(load("./outputs/00_preprocess/list_snrnaseq_sampleid.rda"))[["validation"]]
validation.sample.id <- unique(gsub("_.*", "", sample.id.validate))

# load mae, subset validation
mae.path <- file.path("outputs", "01_mae", "mae_analysis_append.rda")
mae <- get(load(mae.path))
# filter validation ids
filter.mae <- colData(mae)$sample.id %in% sample.id.validate
mae <- mae[,filter.mae,]
unique(colData(mae)[,1][complete.cases(mae)])

# get bulk expression
rse.counts <- mae[["bulk.rnaseq"]]
names(assays(rse.counts)) <- "counts"
# snrnaseq reference -- using same reference across experiments
sce.iter <- mae[["snrnaseq.k2.all"]]
sce.iter <- sce.iter[,sce.iter$Sample %in% sample.id.validate]
sce.iter <- logNormCounts(sce.iter)

# define the sample id vector
sample.id.vector <- unique(sce.iter$Sample)

# define the true proportions
list.df.true <- metadata(sce.iter)[["list.df.true.k2"]]

#--------------
# prep s-vector
#--------------
s.manual <- c("glial" = 3, "neuron" = 10)
s.null <- c("glial" = 1, "neuron" = 1)
list.s.null <- lapply(sample.id.vector, function(id){s.null})
list.s.manual <- lapply(sample.id.vector, function(id){s.manual})
names(list.s.null) <- names(list.s.manual) <- sample.id.vector
list.s.pred <- list(s.null = list.s.null, s.manual = list.s.manual)

#-------------
# run a/b test
#-------------
# experiment variables
deconvolution.algorithm <- "nnls"
cell.type.variable <- "k2"
base.path <- file.path("scripts/05_bulk/abtest/")

# prep dependencies
source(file.path(base.path, "00_helper-functions.R"))
# run experiments
source(file.path(base.path, "02-01-02_rse-counts_logcounts-lutearg_shared-reference-experiments.R"))
source(file.path(base.path, "02-02-02_rse-counts_lognorm-yz_within-reference-experiments.R"))

# prepare results df
shared.colnames <- intersect(colnames(df.s.k2.shared.counts.logcounts),
                             colnames(df.s.k2.within.counts.lognorm))
list.logcounts <- list(
  df.s.k2.shared.counts.logcounts[,shared.colnames],
  df.s.k2.within.counts.lognorm
)
cd <- colData(mae[["bulk.rnaseq"]])
list.logcounts <- lapply(list.logcounts, function(df){
  df$bulk.sample.id <- rownames(df)
  df$bulk.sample.condition <- cd[df$bulk.sample.id,]$expt_condition
  df$assay.name.lutearg <- "logcounts"
  df
})
df.k2 <- as.data.frame(do.call(rbind, list.logcounts))

# append true values
sample.id.vector <- unique(df.k2$sample.id)
df.k2$true.neuron <- df.k2$true.glial <- "NA"
for(sample.id in sample.id.vector){
  df.k2[df.k2$sample.id==sample.id,]$true.neuron <-
    as.numeric(list.df.true[[sample.id]][["neuron"]])
  df.k2[df.k2$sample.id==sample.id,]$true.glial <-
    as.numeric(list.df.true[[sample.id]][["glial"]])
}
df.k2$true.neuron <- as.numeric(df.k2$true.neuron)
df.k2$abs.error.neuron <- abs(as.numeric(df.k2$neuron)-as.numeric(df.k2$true.neuron))
df.k2$abs.error.glial <- abs(as.numeric(df.k2$glial)-as.numeric(df.k2$true.glial))

# clear cache
rm(df.s.k2.shared.counts.logcounts)
rm(df.s.k2.within.counts.lognorm)
rm(mae)
rm(mae.final)
rm(rse.counts)
rm(rse.rpkm)
rm(sce.iter)
# append crossvalidation info
df.k2.validate <- df.k2
df.k2.train <- get(load("./outputs/05_bulk/k2_results.rda"))
filter.tain <- df.k2.train$bulk.scale.type=="counts"
filter.train <- df.k2.train$assay.name.lutearg=="logcounts"
df.k2.train <- df.k2.train[filter.train,]
dim(df.k2.train)
df.k2.train$crossvalidation <- "train"
df.k2.validate$crossvalidation <- "validation"
identical(colnames(df.k2.validate), colnames(df.k2.train))
df.k2 <- as.data.frame(rbind(df.k2.validate, df.k2.train))
dim(df.k2)
rm(df.k2.train)
rm(df.k2.validate)

#-----
# save
#-----

# save
save(df.k2, file = "./outputs/05_bulk/validation_results.rda")

# save environment
save.image(file = "./env/05_bulk/02_crossvalidate_script.RData")