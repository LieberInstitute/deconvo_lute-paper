#!/usr/bin/env R

# Author: Sean Maden
#
# RMSE and nuclei proportions from pseudobulks.
#
#

#------
# load
#------
load("env/08_improvement/02_simulation_script.RData")
df.k2 <- dfPseudobulk
colnames(df.k2) <- c("non.pred", "plasma.pred", 
                     "condition", "plasma.true", "non.true", "sample.id")

#------
# rmse
#------
rmseType <- function(df, type){
  errorVector <- unlist(lapply(type, function(ti){
    catchTrue <- grepl(ti,colnames(df))&grepl("\\.true$",colnames(df))
    catchPred <- grepl(ti,colnames(df))&grepl("\\.pred$",colnames(df))
    df[,catchTrue]-df[,catchPred]
  })) |> as.numeric()
  return(
    sqrt(mean(errorVector^2))
  )
}

rmsePlasmaK2 <- rmseType(df.k2, "plasma") # 0.03796989
rmseNonK2 <- rmseType(df.k2, "non") # 0.03796989
rmseAllK2 <- rmseType(df.k2, c("plasma", "non")) # 0.03796989

#-------------
# proportions
#-------------
summary(df.k2$plasma.true)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000900 0.001700 0.002200 0.003367 0.005175 0.008000
summary(df.k2$non.true)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.9920  0.9948  0.9978  0.9966  0.9983  0.9991

#-------------------
# nuclei per sample
#-------------------
mrna.yield <- 
  read.table("./data/monaco_et_al_2019/shared/mRNAyield_PBMC.txt", header = T)
as.numeric(mrna.yield$Cell_no.) |> 
  median() |> format(scientific = T)

#----------------------
# fold size difference
#----------------------
as.numeric(cellSizes["Plasmablasts"])/
  as.numeric(mean(cellSizes[!names(cellSizes)=="Plasmablasts"]))
# 15.31991



