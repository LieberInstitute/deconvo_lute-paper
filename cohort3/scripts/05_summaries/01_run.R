#!/usr/bin/env R

# Author: Sean Maden
#
# Summarize dataset elements.
# 
# * "GSE107011_Processed_data_TPM.txt", GEO record.
#
#
#
#

libv <- c("SummarizedExperiment", "biomaRt")
sapply(libv, library, character.only = TRUE)
source("./source/00_read_experiment_data.R")
experimentData <- getExperimentData()

#--------------------------------------------
# flow cytometry "true" proportions summaries
#--------------------------------------------
# plasmablasts
proportionsPlasmablast <- 
  experimentData$p.true["Plasmablasts",] |> as.numeric() |> summary()
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0007014 0.0014741 0.0020032 0.0031640 0.0049517 0.0077615

# non-plasmablasts
proportionsNonplasmablastsCombined <- 
  summary(1-as.numeric(experimentData$p.true["Plasmablasts",]))
proportionsNonplasmablastsCombined
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.9922  0.9950  0.9980  0.9968  0.9985  0.9993

# range of non-plasmablast cell composition
proportionsNonplasmablastsAlltypes <- experimentData$p.true
proportionsNonplasmablastsAlltypes <-
  pTrue[!rownames(pTrue)=="Plasmablasts",] |> unlist() |> as.numeric() |> summary()
proportionsNonplasmablastsAlltypes
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.001049 0.015634 0.036024 0.062302 0.088761 0.275502

proportionsTable <- rbind(proportionsPlasmablast, 
                          rbind(proportionsNonplasmablastsCombined, 
                                proportionsNonplasmablastsAlltypes)) |>
  as.data.frame()
rownames(proportionsTable) <- 
  c("plasmablast", "nonplasmablast", "nonplasmablast_unbinned")

#-----
# save
#-----
# output proportions table
save(proportionsTable, 
     file = "./outputs/08_improvements/cellProportionsTable.rda")
data.table::fwrite(
  proportionsTable,"./outputs/08_improvements/cellProportionsTable.csv")

save.image("./env/05_summaries/01_run_script.RData")

