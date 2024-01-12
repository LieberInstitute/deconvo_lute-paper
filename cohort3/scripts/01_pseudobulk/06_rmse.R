#!/usr/bin/env R

# Author: Sean Maden
#
# RMSE and nuclei proportions from pseudobulks. Uses rmseTable().
#
#

#------
# load
#------
load("env/08_improvement/02_simulation_script.RData")
df.k2 <- dfPseudobulk
#colnames(df.k2) <- c("non.pred", "plasma.pred", 
#                     "condition", "plasma.true", "non.true", "sample.id")

colnames(df.k2) <- c("NonPlasmablasts.pred", "Plasmablasts.pred", 
                     "condition", "Plasmablasts.true", "NonPlasmablasts.true",
                     "sample.id")

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

rmseTableK2 <- function(df.k2, cohortLabel){
  # rmseTable
  # df.k2 table for k2 results
  # cohortLabel dataset cohort label
  #
  # returns rmse table for k2 groups.
  #
  
  #-------------------
  # calculations -- K2
  #-------------------
  # K2
  # SCALING
  listRmseK2 <- list(
    c("NonPlasmablasts"), c("Plasmablasts"), 
    c("NonPlasmablasts", "Plasmablasts")
  )
  rmseK2Scale <- lapply(listRmseK2, function(typeVector){
    rmseType(df.k2[df.k2$condition=="scaled",], typeVector)
  })
  rmseK2Scale
  
  names(rmseK2Scale) <- sapply(listRmseK2, function(i){paste0(i,collapse =";")})
  
  # NO SCALING
  # SCALING
  rmseK2Noscale <- lapply(listRmseK2, function(typeVector){
    round(rmseType(df.k2[df.k2$condition=="unscaled",], typeVector),3)
  })
  names(rmseK2Noscale) <- sapply(listRmseK2, function(i){paste0(i,collapse =";")})
  
  #-----------
  # SUPPLEMENT
  #-----------
  
  returnListScale <- list(c(rmseK2Scale))
  returnTableScale <- do.call(cbind, lapply(seq(length(returnListScale)), function(ii){
    c(names(returnListScale)[ii], returnListScale[[ii]])
  })) |> as.data.frame()
  returnTableScale$k <- c(rep("k2", length(listRmseK2)))
  
  returnListNoscale <- list(c(rmseK2Noscale))
  returnTableNoscale <- do.call(cbind, lapply(seq(length(returnListNoscale)), function(ii){
    c(names(returnListNoscale)[ii], returnListNoscale[[ii]])
  })) |> as.data.frame()
  returnTableNoscale$k <- c(rep("k2", length(listRmseK2)))
  
  returnTableScale$type <- "scale"
  returnTableNoscale$type <- "noscale"
  returnTable <- rbind(returnTableScale, returnTableNoscale)
  colnames(returnTable)[1] <- "rmse"
  
  returnTable$cohort <- cohortLabel
  returnTable$numberCellTypes <- sapply(
    rownames(returnTable), function(ii){
      length(unlist(strsplit(ii, "\\;")))}) |> as.numeric()
  
  return(returnTable)
}

#----------------
# make rmse table
#----------------
supplementTable <- rmseTableK2(df.k2, "cohort3")

#-----
# save
#-----
# save
save(supplementTable, file = "./outputs/08_improvements/rmse_supplementTable.rda")


