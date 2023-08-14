#!/usr/bin/env R

# Author: Sean Maden
#
#
#

libv <- c("MultiAssayExperiment")
sapply(libv, library, character.only = TRUE)

# load
mae.name <- "mae_with-rpkm_additional-data_final.rda"
mae.path <- file.path("deconvo_method-paper", "outputs", "01_prepare-datasets", "mae.name")
mae <- get(load(mae.path))

# get assay labels
labels <- names(mae)
labels
#[1] "sn1.rnaseq"           "sn2.rnaseq"           "sn3.rnaseq"           "bulk.rnaseq"          "bulk.rpkm.rnaseq"    
#[6] "rnascope.image"       "df.cellstat.rnascope"

# get training sample ids from sce marker data
sce.train.sample.id <- unique(gsub("_.*" ,"",mae[["sn1.rnaseq"]]$Sample))
# sce.train.sample.id
# [1] "Br2720" "Br6471" "Br8492" "Br2743" "Br3942" "Br6423" "Br8325"
length(sce.train.sample.id)
# [1] 7

# get training ybulk sample ids
ybulk.train.sample.id.all <- unique(gsub("_.*" ,"",mae[["bulk.rnaseq"]]$BrNum))
ybulk.train.sample.id <- ybulk.train.sample.id.all[ybulk.train.sample.id.all %in% sce.train.sample.id]
ybulk.train.sample.id
# [1] "Br2720" "Br6471" "Br2743" "Br3942" "Br6423" "Br8325" "Br8492"
length(ybulk.train.sample.id)
# [1] 7

# get test/validate ybulk sample ids
ybulk.test.validate.sample.id <- ybulk.train.sample.id.all[!ybulk.train.sample.id.all %in% sce.train.sample.id]
# ybulk.test.validate.sample.id
# [1] "Br6432" "Br6522" "Br8667"

# save mae validate subset
# test/validate subset
train.string <- paste0(ybulk.train.sample.id.all, collapse = "|")
train.mae.filter <- grepl(train.string, mae$sample.id)
mae.validate <- mae[,train.mae.filter,]
# save
mae.validate.name <- "mae-validate_with-rpkm_cohort1.rda"
mae.validate.path <- file.path("deconvo_method-paper", "outputs", "01_prepare-datasets", "mae.name")
save(mae.validate, file = mae.validate.path)






