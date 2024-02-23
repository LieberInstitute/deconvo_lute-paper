#!/usr/bin/env R

# Author: Sean Maden
#
# Summarize a dataset
#

get_data_summaries <- function(df, var.plot){
  list(dim = dim(df),
       summary = summary(df[,var.plot]),
       unique.varplot = unique(df[,var.plot]))
}