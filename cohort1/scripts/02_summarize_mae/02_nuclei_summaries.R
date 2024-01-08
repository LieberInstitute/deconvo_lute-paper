#!/usr/bin/env R

# Author: Sean Maden
#
# Summarizes nuclei, cohort1.
#
#

#source("./cohort2/scripts/01_pseudobulk/00_parameters-pseudobulk.R")
#sapply(libv, library, character.only = T)

#------
# load
#------
# load marker data
load("env/02_pseudobulk/01_k2.RData")

sceK2 <- list.sce.markers[["k2"]]
sceK3 <- list.sce.markers[["k3"]]

dfSummary <- function(sce, 
                      sampleIdVariable = "donor",
                      variablesVector = c("k2", "k3"),
                      typeNames = c("neuron", "glial", "Inhib", "Excit")){
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
        summary(markerSizes)
        
        c(variable,
          type,
          median(sampleVector), # cells/nuclei
          mean(sampleVector), # cells/nuclei
          sd(sampleVector), # cells/nuclei
          median(proportionsVector), # proportion per sample
          mean(proportionsVector), # proportion per sample
          sd(proportionsVector), # proportion per sample
          mean(markerSizes), # marker library size
          median(markerSizes), # marker library size
          sd(markerSizes) # marker library size
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

summaryK2 <- dfSummary(sceK2, "Sample", "k2", c("neuron", "glial"))
summaryK3 <- dfSummary(sceK3, "Sample", "k3", c("Inhib", "Excit", "glial"))
summaryTable <- rbind(summaryK2, summaryK3)


#-----
# save
#-----
data.table::fwrite(summaryTable, file = "./outputs/02_summarize_mae/cohort1_nuclei_summaries.csv")
