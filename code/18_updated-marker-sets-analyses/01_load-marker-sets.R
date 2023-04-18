# load the main marker sets from using 2 approaches in dlpfc ro1 training

libv <- c("here", "dplyr", "lute", "biomaRt", "SummarizedExperiment", 
          "SingleCellExperiment")
sapply(libv, library, character.only = TRUE)

# load k2 markers expression
outputs.path <- here("deconvo_method-paper", "outputs", 
                     "18_updated-marker-sets-analyses")
# ro1 train
se.train.path <- here(outputs.path, "se-adj_marker-filter-all_train.rda")
se.train <- get(load(se.train.path))
marker.train.overlaps.path <- here(outputs.path, "overlap-markers-filtered-k2_train.rda")
marker.train.overlaps <- get(load(marker.train.overlaps.path))
marker.train.group.path <- here(outputs.path, "group-markers-k2_train.rda")
marker.train.group <- get(load(marker.train.group.path))
# mrb
se.mrb.path <- here(outputs.path, "se-adj_marker-filter-all_mrb.rda")
se.mrb <- get(load(se.mrb.path))
marker.mrb.overlaps.path <- here(outputs.path, "overlap-markers-filtered-k2_mrb.rda")
marker.mrb.overlaps <- get(load(marker.train.overlaps.path))

# sutton markers
format.name <- "sutton-brain-markers-formatted.rda"
format.path <- here("deconvo_method-paper", "outputs", 
                    "18_updated-marker-sets-analyses", format.name)
sutton <- get(load(format.path))
