#!/usr/bin/env R

# Author: Sean Maden
#
# RMSE from pseudobulks, cohort2.
#
#

# load
load("./env/01_pseudobulk/01_k2_mrb_script.RData")
df.k2 <- dfp.tall
rm(dfp.tall)

load("env/01_pseudobulk/01_k3_mrb_script.RData")
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

#----------------
# make rmse table
#----------------
rmseTable <- function(df.k2, df.k3, cohortLabel){
  # rmseTable
  # df.k2 table for k2 results
  # df.k3 table for k3 results
  # cohortLabel dataset cohort label
  #
  # returns rmse table for k2 and k3 groups.
  #
  
  #-------------------
  # calculations -- K2
  #-------------------
  # K2
  # SCALING
  listRmseK2 <- list(
    c("neuron"), c("glial"), c("neuron", "glial")
  )
  rmseK2Scale <- lapply(listRmseK2, function(typeVector){
    rmseType(df.k2[df.k2$type=="withscale",], typeVector)
  })
  names(rmseK2Scale) <- sapply(listRmseK2, function(i){paste0(i,collapse =";")})
  
  # NO SCALING
  # SCALING
  rmseK2Noscale <- lapply(listRmseK2, function(typeVector){
    round(rmseType(df.k2[df.k2$type=="noscale",], typeVector),3)
  })
  names(rmseK2Noscale) <- sapply(listRmseK2, function(i){paste0(i,collapse =";")})
  
  #-------------------
  # calculations -- K3
  #-------------------
  # WITH SCALING
  # NO SCALE
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
  
  # WITH SCALE
  rmseK3Scale <- lapply(listRmseK3, function(typeVector){
    round(rmseType(df.k3[df.k3$type=="withscale",], typeVector), 3)
  })
  names(rmseK3Scale) <- sapply(listRmseK3, function(i){paste0(i,collapse =";")})
  # NO SCALE
  rmseK3Noscale <- lapply(listRmseK3, function(typeVector){
    round(rmseType(df.k3[df.k3$type=="noscale",], typeVector), 3)
  })
  names(rmseK3Noscale) <- sapply(listRmseK3, function(i){paste0(i,collapse =";")})
  
  
  #-----------
  # SUPPLEMENT
  #-----------
  
  returnListScale <- list(c(rmseK2Scale, rmseK3Scale))
  returnTableScale <- do.call(cbind, lapply(seq(length(returnListScale)), function(ii){
    c(names(returnListScale)[ii], returnListScale[[ii]])
  })) |> as.data.frame()
  returnTableScale$k <- c(rep("k2", length(listRmseK2)),rep("k3", length(listRmseK3)))
  
  returnListNoscale <- list(c(rmseK2Noscale, rmseK3Noscale))
  returnTableNoscale <- do.call(cbind, lapply(seq(length(returnListNoscale)), function(ii){
    c(names(returnListNoscale)[ii], returnListNoscale[[ii]])
  })) |> as.data.frame()
  returnTableNoscale$k <- c(rep("k2", length(listRmseK2)),rep("k3", length(listRmseK3)))
  
  returnTableScale$type <- "scale"
  returnTableNoscale$type <- "noscale"
  returnTable <- rbind(returnTableScale, returnTableNoscale)
  colnames(returnTable)[1] <- "rmse"
  
  returnTable$cohort <- cohortLabel
  returnTable$numberCellTypes <- sapply(
    rownames(returnTable), function(ii){
      length(unlist(strsplit(ii, "\\.")))}) |> as.numeric()
  
  return(returnTable)
}

supplementTable <- rmseTable(df.k2, df.k3, "cohort2")

# gets small k3 values
rmseType(df.k3[df.k3$type=="withscale",], "glial") # 1.433292e-16
rmseType(df.k3[df.k3$type=="withscale",], "Inhib") # 1.488227e-16
rmseType(df.k3[df.k3$type=="withscale",], "Excit") # 2.028564e-16
rmseType(df.k3[df.k3$type=="withscale",], "neuron") # 1.155557e-16
rmseType(df.k3[df.k3$type=="withscale",], c("glial","Excit")) # 1.756331e-16
rmseType(df.k3[df.k3$type=="withscale",], c("glial","Inhib")) # 1.461018e-16
rmseType(df.k3[df.k3$type=="withscale",], c("Excit","Inhib")) # 1.779029e-16
rmseType(df.k3[df.k3$type=="withscale",], c("Excit","Inhib", "glial")) # 1.671747e-16

#-----
# save
#-----
# save
save(supplementTable, file = "./outputs/01_pseudobulk/rmse_supplementTable.rda")