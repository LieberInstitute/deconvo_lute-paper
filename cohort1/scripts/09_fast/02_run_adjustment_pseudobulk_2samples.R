source("./scripts/08_adjustment/00_sopt.R")
source("./scripts/08_adjustment/00_sopt_utilities.R")
source("./scripts/08_adjustment/00_param.R")
libv <- c("snow", "dplyr", "parallel", "doParallel", "lute", "BisqueRNA", "MuSiC", 
          "dplyr", "MultiAssayExperiment", "GGally")
sapply(libv, library, character.only = T)
source("./scripts/08_adjustment/00_musicParam-class.R")
num.dfs.steps <- 5
mae <- get(load("outputs/01_mae/mae_analysis_append.rda"))
sample.id.keep <- c("Br8325_mid", "Br3942_mid")
mae <- mae[,colData(mae)$sample.id %in% sample.id.keep,]

#mae[["bulk.pb.k2"]] <- cbind(mae[["bulk.pb.k2"]],mae[["bulk.pb.k2"]])
list.experiment.results <- experiment_all_samples(
  colData(mae)$sample.id, 
  mae, 
  bulk.mae.name = "bulk.pb.k2", 
  dfs.steps = num.dfs.steps)

df.res <- as.data.frame(
  do.call(rbind, lapply(list.experiment.results, function(item){item$df.res})))
df.res$sample.id <- gsub("_.*", "", rownames(df.res))
list.dfp <- get_dfp_list(df.res)

# save image
rm(mae)
save.image(file = "./env/09_fast/02_run_2sample_script.RData")
