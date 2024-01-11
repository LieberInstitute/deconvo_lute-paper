#!/usr/bin/env R

# Author: Sean Maden
#
# Summarizes samples, cohort2.
#
#

libv <- c("SingleCellExperiment")
sapply(libv, library, character.only = TRUE)

# load
load("./env/01_pseudobulk/01_k2_mrb_script.RData")

#--------
# summary
#--------
cd <- colData(sce)

# nuclei by donor
dfs <- table(cd$donor) |> as.data.frame()
summary(dfs[,2])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1686    2948    4209    3722    4740    5270

# proportions nuclei by k2 type
dfs <- table(cd$donor, cd$k2)
dfp <- apply(dfs, 2, function(ci){ci/sum(ci)})
summary(dfp[,"neuron"])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1855  0.1997  0.2140  0.3333  0.4073  0.6006

summary(dfp[,"glial"])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1163  0.2587  0.4011  0.3333  0.4419  0.4826

sampleSummaries <- function(sce, round = 3){
  #
  # sce SingleCellExperiment
  #
  
  cd <- colData(sce)
  
  round(, digits = 3)
  
  summariesTable <- c(
    
  )
  return(summariesTable)
}


#-----
# save
#-----

