# load the main marker sets from using 2 approaches in dlpfc ro1 training

libv <- c("here", "dplyr", "lute", "biomaRt", "SummarizedExperiment", 
          "SingleCellExperiment")
sapply(libv, library, character.only = TRUE)

# load markers -- training
labels.vector <- c("k2", "k3", "k4")
# get approach1 markers
list.name <- "list-scef_markers-k2-k3-k4_ro1-dlpfc.rda"
list.path <- here("deconvo_method-paper", "outputs", "09_manuscript", list.name)
sce.list1 <- get(load(list.path))
marker.list1 <- lapply(sce.list1, rownames)
names(marker.list1) <- labels.vector

# get approach2 markers
load.name.vector <- paste0("overall-markers-filtered-",labels.vector,"_train.rda")
load.path <- here("deconvo_method-paper", "outputs", 
                  "15_k2-simulations_within-sample-matched", load.name.vector)
marker.list2 <- lapply(load.path, function(path){get(load(path))})
names(marker.list2) <- labels.vector

# load markers -- mrb
# approach1 
list.name <- "list-scef_markers-k2-k3-k4_mrb-dlpfc.rda"
list.path <- here("deconvo_method-paper", "outputs", "09_manuscript", list.name)
sce.mrb.list1 <- get(load(list.path))
marker.mrb.list1 <- lapply(sce.mrb.list1, rownames)
names(marker.mrb.list1) <- labels.vector
# approach2
k2.name <- "k2-markers_overall-filtered_mrb.rda"
k2.path <- here("deconvo_method-paper", "outputs", 
                "15_k2-simulations_within-sample-matched", k2.name)
k2.mrb <- get(load(k2.path))
# load sce
# do marker selection

# sutton markers
format.name <- "sutton-brain-markers-formatted.rda"
format.path <- here("deconvo_method-paper", "outputs", 
                    "18_updated-marker-sets-analyses", format.name)
sutton <- get(load(format.path))
