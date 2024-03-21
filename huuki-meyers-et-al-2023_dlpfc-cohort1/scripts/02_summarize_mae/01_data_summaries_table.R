#!/usr/bin/env R

# Author: Sean Maden
#
# Make data summaries table for manuscript.
#
#
#
#
#

library(dplyr)

source("./scripts/02_summarize_mae/00_param.R")

# load
# dataset names
path.data.pre.filter <- "./outputs/01_mae/mae_allsamples_append.rda"
path.data.post.filter <- "./outputs/01_mae/mae_analysis_append.rda"
mae.all <- get(load(path.data.pre.filter))
mae.filter <- get(load(path.data.post.filter))

list.sample.id <- get(load("./outputs/00_preprocess/list_snrnaseq_sampleid.rda"))
sample.id.train <- list.sample.id[["train"]]
sample.id.validate <- list.sample.id[["validation"]]

# summaries
list.summaries.prefilter <- summaries_df_list(mae.all, filter.type.label = "prefilter")
list.summaries.postfilter <- summaries_df_list(mae.filter, filter.type.label = "postfilter")

# write summary tables
write.csv(list.summaries.prefilter[["wide"]], 
          file = "./outputs/02_summarize_mae/TABLE_S8_SUBSET.csv", row.names = FALSE)

# save image
save.image(file = './env/02_summarize_mae/01_data_summaries_script.RData')
