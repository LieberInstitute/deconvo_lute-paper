#!/usr/bin/env R

# Author: Sean Maden
#
# Run successive iterations of S optimization. Use the real bulk data.
#

libv <- c("MultiAssayExperiment", "ggplot2", "gridExtra")
sapply(libv, library, character.only = T)

source("./source/00_dataset_summaries.R")
source("./source/00_deconvo_plots.R")
source("./source/00_sopt.R")
# source("./source/00_sopt.R")

# load
mae <- get(load("./outputs/01_mae/mae_analysis_append.rda"))
mae <- mae[,colData(mae)$sample.id=="Br8492_mid",]

# new run
dfs.param <- data.frame(
  s.min = c(1, NA, NA, NA),
  s.max = c(80, NA, NA, NA),
  s.step = c(2, 1e-2, 1e-3, 1e-4),
  s.diff = c(NA, 1e-1, 1e-2, 1e-3)
)
sample.id.vector <- "Br8492_mid"
df.true.list <- metadata(mae[[1]])[["list.df.true.k2"]]
y.unadj <- mae[["bulk.rnaseq"]][,1]
sce <- mae[[1]]
list.res <- run_sopt_series(dfs.param, sample.id.vector = sample.id.vector, 
                            df.true.list = df.true.list, y.unadj = y.unadj, 
                            sce = sce)

identical(list.res$iter1$df.res,
          list.res$iter2$df.res)

list.res.plots <- plot_sopt_series(list.res)

save.image(file = "./env/09_fast/01_run_bulk_1sample_script.RData")
