#!/usr/bin/env R

# Author: Sean Maden
#
# RMSE from pseudobulks, cohort1.
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

#--------------------
# append neuron to k3
#--------------------
df.k3$neuron.true <- NA
for(sample.id in unique(df.k2$sample.id)){
  if(sample.id %in% df.k3$sample.id){
    df.k3[df.k3$sample.id==sample.id,]$neuron.true <- 
      df.k2[df.k2$sample.id==sample.id,]$neuron.true[1]
  }
}

# append neuron pred as sum
df.k3$neuron.pred <- df.k3$Excit.pred + df.k3$Inhib.pred

#-------------
# calculations
#-------------
listRmseK2 <- list(
  c("neuron"), c("glial"), c("neuron", "glial")
)
rmseK2 <- lapply(listRmseK2, function(typeVector){
  round(rmseType(df.k2, typeVector),3)
})
names(rmseK2) <- sapply(listRmseK2, function(i){paste0(i,collapse =";")})
rmseK2

listRmseK3 <- list(
  c("neuron"), 
  c("glial"), 
  c("Excit"), 
  c("Inhib"), 
  c("Excit", "glial"), 
  c("Inhib", "glial"), 
  c("Excit", "Inhib"), 
  c("Excit", "Inhib", "glial")
)
rmseK3 <- lapply(listRmseK3, function(typeVector){
  round(rmseType(df.k3, typeVector), 3)
})
names(rmseK3) <- sapply(listRmseK3, function(i){paste0(i,collapse =";")})
rmseK3