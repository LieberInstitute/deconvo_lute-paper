#!/usr/bin/env R

# Author: Sean Maden
#
# RMSE from pseudobulks.
#
#

# load
load("env/02_pseudobulk/01_k2.RData")
df.k2 <- dfp.tall
rm(dfp.tall)

load("env/02_pseudobulk/02_k3.RData")
df.k3 <- dfp.tall
rm(dfp.tall)

# rmse
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

rmseNeuronK2 <- rmseType(df.k2, "neuron")
rmseGlialK2 <- rmseType(df.k2, "glial")
rmseAllK2 <- rmseType(df.k2, c("neuron", "glial")) # 0.15
rmseAllK3 <- rmseType(df.k3, c("Excit", "Inhib", "glial")) # 0.11
