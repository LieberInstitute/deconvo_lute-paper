#!/usr/bin/env R

# Author: Sean Maden
#
# Summarizes nuclei, cohort2.
#
#

#source("./cohort2/scripts/01_pseudobulk/00_parameters-pseudobulk.R")
#sapply(libv, library, character.only = T)

#------
# load
#------
# load marker data
load("env/01_pseudobulk/01_k2_mrb_script.RData")
load("env/01_pseudobulk/01_k3_mrb_script.RData")

dfSummary <- function(sce, 
                      sampleIdVariable = "donor",
                      variablesVector = c("k2", "k3"),
                      typeNames = c("neuron", "glial", "Inhib", "Excit"),
                      roundDigits = 3){
  # dfSummary
  #
  # sce SingleCellExperiment
  # variablesVector vector of variables
  # typeNames vector of celltypes
  #
  cd <- colData(sce)
  listSummaryValues <- lapply(variablesVector, function(variable){
    # proportions
    donorTableTotal <- table(cd[,sampleIdVariable]) |> as.data.frame()
    variableValuesList <- lapply(typeNames, function(type){
      message(type, "; ", variable)
      
      if(type %in% cd[,variable]){
        cdf <- cd[cd[,variable]==type,]
        sampleCellCountsVector <- 
          table(cdf[,variable], cdf[,sampleIdVariable]) |> as.data.frame()
        sampleVector <- sampleCellCountsVector[,3]
        proportionsVector <- 
          sampleVector/donorTableTotal[,2]
        
        markerSizes <- sce[,colData(sce)[,variable]==type] |> 
          counts() |> colSums()
        
        c(variable,
          type,
          round(median(sampleVector), roundDigits), # cells/nuclei
          round(mean(sampleVector), roundDigits), # cells/nuclei
          round(sd(sampleVector), roundDigits), # cells/nuclei
          round(median(proportionsVector), roundDigits), # proportion per sample
          round(mean(proportionsVector), roundDigits), # proportion per sample
          round(sd(proportionsVector), roundDigits), # proportion per sample
          round(mean(markerSizes), roundDigits), # marker library size
          round(median(markerSizes), roundDigits), # marker library size
          round(sd(markerSizes), roundDigits) # marker library size
        )
      }
      
    })
    do.call(rbind, variableValuesList) |> as.data.frame()
  })
  dfSummaries <- do.call(rbind, listSummaryValues) |> as.data.frame()
  colnames(dfSummaries) <- c("k", "cellType", 
                             "median_nuclei_per_sample",
                             "mean_nuclei_per_sample",
                             "sd_nuclei_per_sample",
                             "median_proportion_nuclei_per_sample",
                             "mean_proportion_nuclei_per_sample",
                             "sd_proportion_nuclei_per_sample",
                             "median_marker_library_size",
                             "mean_marker_library_size",
                             "sd_marker_library_size")
  return(dfSummaries)
}

summaryK2 <- dfSummary(sce.k2, "donor", "k2", c("neuron", "glial"))
summaryK3 <- dfSummary(sce.k3, "donor", "k3", c("Inhib", "Excit", "glial"))
summaryTable <- rbind(summaryK2, summaryK3)

#-----
# save
#-----
data.table::fwrite(summaryTable, file = "./outputs/02_summaries/cohort2_nuclei_summaries.csv")

