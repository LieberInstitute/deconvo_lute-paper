# load
supplementTable1 <- 
  get(load("./outputs/01_pseudobulk/rmse_supplementTable.rda"))
supplementTable2 <- 
  get(load("./outputs/02_pseudobulk/rmse_supplementTable.rda"))
# bind

tableOut <- 
  rbind(supplementTable1, supplementTable2, supplementTable3) |> as.data.frame()

# save
write.csv("./out/03_rmse/rmse_supplementTable.csv")