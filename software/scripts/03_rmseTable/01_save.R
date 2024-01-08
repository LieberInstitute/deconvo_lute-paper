# load
# setwd("..") # set path to deconvo_lute-paper
supplementTable1 <- 
  get(load("./cohort1/outputs/02_pseudobulk/rmse_supplementTable.rda"))
supplementTable2 <- 
  get(load("./cohort2/outputs/01_pseudobulk/rmse_supplementTable.rda"))
supplementTable3 <- 
  get(load("./cohort3/outputs/08_improvements/rmse_supplementTable.rda"))

# bind
tableOut <- 
  rbind(supplementTable1, supplementTable2, supplementTable3) |> as.data.frame()

tableOut$cellTypes <- rownames(tableOut)

# save
# write.csv(tableOut, file = "./software/outputs/03_rmseTable/rmse_supplementTable.csv")
data.table::fwrite(tableOut, file = "./software/outputs/03_rmseTable/rmse_supplementTable.csv")