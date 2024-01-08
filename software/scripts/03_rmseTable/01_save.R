# load
# setwd("..") # set path to deconvo_lute-paper
supplementTable1 <- 
  get(load("./cohort1/outputs/02_pseudobulk/rmse_supplementTable.rda"))
supplementTable2 <- 
  get(load("./cohort2/outputs/01_pseudobulk/rmse_supplementTable.rda"))
supplementTable3 <- 
  get(load("./cohort2/outputs/01_pseudobulk/rmse_supplementTable.rda"))

# bind

tableOut <- 
  rbind(supplementTable1, supplementTable2, supplementTable3) |> as.data.frame()

# save
write.csv("./software/out/03_rmse/rmse_supplementTable.csv")