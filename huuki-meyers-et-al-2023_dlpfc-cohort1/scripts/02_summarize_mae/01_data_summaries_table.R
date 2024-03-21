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

# save image
save.image(file = './env/02_summarize_mae/01_data_summaries_script.RData')
