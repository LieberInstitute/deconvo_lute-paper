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

# dataset names
path.data.pre.filter <- "./outputs/01_mae/mae_allsamples_append.rda"
path.data.post.filter <- "./outputs/01_mae/mae_analysis_append.rda"

list.sample.id <- get(load("./outputs/00_preprocess/list_snrnaseq_sampleid.rda"))
sample.id.train <- list.sample.id[["train"]]
sample.id.validate <- list.sample.id[["validation"]]

source("./scripts/02_summarize_mae/00_param.R")

# load
mae.all <- get(load(path.data.pre.filter))
mae.filter <- get(load(path.data.post.filter))

# summaries
list.summaries.prefilter <- summaries_df_list(mae.all, filter.type.label = "prefilter")
list.summaries.postfilter <- summaries_df_list(mae.filter, filter.type.label = "postfilter")

# write summary tables
write.csv(list.summaries.prefilter[["df.wide"]], 
          file = "./outputs/02_summarize_mae/table2_platforms.csv", row.names = FALSE)
df.filter.all <- rbind(list.summaries.postfilter[[1]],
                       list.summaries.prefilter[[1]])
write.csv(df.filter.all, 
          file = "./outputs/02_summarize_mae/table_compare_filters.csv", row.names = FALSE)

# save image
save.image(file = './env/02_summarize_mae/01_data_summaries_script.RData')
