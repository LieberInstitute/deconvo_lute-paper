#!/usr/bin/env R

#
# Get cell type proportions
#




libv <- c("dplyr", "ggplot2")
sapply(libv, library, character.only = TRUE)

folder.name <- "08_sizes"
source(file.path("./source/00_cell_sizes.R"))

# load
mae <- get(load("./outputs/01_mae/mae_analysis_append.rda"))
sn <- mae[["snrnaseq.k2.all"]]
sample.id.variable <- unique(sn[["Sample"]])
mae <- mae[,colData(mae)$sample.id %in% sample.id.variable,]

# sopt cell size factors
sopt <- get(load("./outputs/06_estimate/train_result.rda"))






#----------------
# get cell sizes
#----------------


list.dfp <- list_dfp_wide_tall_size(mae)

dfp.wide <- list.dfp$dfp.wide

dfp.tall <- list.dfp$dfp.tall

# get factors from optimizations
df.sopt <- do.call(rbind, lapply(unique(sopt$sample.id), function(sample.id){
  sf <- sopt[sopt$sample.id==sample.id,]
  c(sample.id, mean(sf$s.glial), mean(sf$s.neuron))
}))

# get centered sizes
dfp.tall$centered.neuron <- scale(dfp.tall$neuron)
dfp.tall$centered.glial <- scale(dfp.tall$glial)

# get ratios
dfp.tall$ratio.glial.neuron <- dfp.tall$glial/dfp.tall$neuron
dfp.tall$ratio.neuron.glial <- dfp.tall$neuron/dfp.tall$glial
dfp.tall$ratio.centered.glial.neuron <- dfp.tall$centered.glial/dfp.tall$centered.neuron
dfp.tall$ratio.centered.neuron.glial <- dfp.tall$centered.neuron/dfp.tall$centered.glial

#------
# save 
#------

# save env
rm(mae)
save.image(file.path("./env/",folder.name,"/01_sizes_script.RData"))
