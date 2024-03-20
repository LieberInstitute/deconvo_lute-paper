#!/usr/bin/env R

# Author: Sean Maden
#
# Build env for rmse TABLE S4. 
# 
# Run from top level (`deconvo_lute-paper/`).
#
#
#

source("software_manuscript/src/rmse.R")

# Parse proportions outputs
# Cohort 1 pseudobulk, K2 and K3
# proportions cohort 1

# k2 results
load(file.path("huuki-meyers-et-al-2023_dlpfc-cohort1/env/02_pseudobulk/01_k2.RData"))
proportionsCohort1K2 <- dfp.tall
rm(dfp.tall)
# k3 results
load(file.path("huuki-meyers-et-al-2023_dlpfc-cohort1/env/02_pseudobulk/02_k3.RData"))
proportionsCohort1K3 <- dfp.tall
rm(dfp.tall)

# k2 results
cellTypesVector <- c("glial", "neuron")
proportionsCohort1K2 <- proportionsCohort1K2[,c(1:6)]
colnames(proportionsCohort1K2) <- c(
  paste0(cellTypesVector, "_known"),
  paste0(cellTypesVector, "_predicted"),
  c("condition.id", "sample.id")
)
# k3 results
cellTypesVector <- c("glial", "Excit", "Inhib")
proportionsCohort1K3 <- proportionsCohort1K3[,c(1:8)]
colnames(proportionsCohort1K3) <- c(
  paste0(cellTypesVector, "_known"),
  paste0(cellTypesVector, "_predicted"),
  c("condition.id", "sample.id")
)

# Cohort 2 pseudobulk, K2 and K3
# proportions cohort 2
# k2 results
load(file.path("./tran-et-al-2021_dlpfc-cohort2/env/01_pseudobulk/01_k2_mrb_script.RData"))
proportionsCohort2K2 <- dfp.tall
rm(dfp.tall)
# k3 results
load(file.path("./tran-et-al-2021_dlpfc-cohort2/env/01_pseudobulk/01_k3_mrb_script.RData"))
proportionsCohort2K3 <- dfp.tall
rm(dfp.tall)

# k2 results
cellTypesVector <- c("glial", "neuron")
proportionsCohort2K2 <- proportionsCohort2K2[,c(1:6)]
colnames(proportionsCohort2K2) <- c(
  paste0(cellTypesVector, "_known"),
  paste0(cellTypesVector, "_predicted"),
  c("condition.id", "sample.id")
)
# k3 results
cellTypesVector <- c("glial", "Excit", "Inhib")
proportionsCohort2K3 <- proportionsCohort2K3[,c(1:8)]
colnames(proportionsCohort2K3) <- c(
  paste0(cellTypesVector, "_known"),
  paste0(cellTypesVector, "_predicted"),
  c("condition.id", "sample.id")
)

# Cohort 1 bulk
load(file.path("./huuki-meyers-et-al-2023_dlpfc-cohort1/env/07_adjustment/03_run_adjustment_realbulk_all_script.RData"))

# nnls
proportionsBulk <-
  data.frame(
    glial_known = 1-dfp.wide$neuron.true,
    neuron_known = dfp.wide$neuron.true,
    glial__predicted = 1-dfp.wide$neuron.nnls.scale,
    neuron_predicted = dfp.wide$neuron.nnls.scale,
    sample.id = rownames(dfp.wide),
    condition.id = rep("scale", nrow(dfp.wide)),
    replicate.id = rep("nnls", nrow(dfp.wide))
  )
proportionsBulk <- rbind(proportionsBulk,
                         data.frame(
                           glial_known = 1-dfp.wide$neuron.true,
                           neuron_known = dfp.wide$neuron.true,
                           glial__predicted = 1-dfp.wide$neuron.nnls.noscale,
                           neuron_predicted = dfp.wide$neuron.nnls.noscale,
                           sample.id = rownames(dfp.wide),
                           condition.id = rep("noscale", nrow(dfp.wide)),
                           replicate.id = rep("nnls", nrow(dfp.wide))
                         ))

# music
proportionsBulk <- rbind(proportionsBulk,
                         data.frame(
                           glial_known = 1-dfp.wide$neuron.true,
                           neuron_known = dfp.wide$neuron.true,
                           glial__predicted = 1-dfp.wide$neuron.music.scale,
                           neuron_predicted = dfp.wide$neuron.music.scale,
                           sample.id = rownames(dfp.wide),
                           condition.id = rep("scale", nrow(dfp.wide)),
                           replicate.id = rep("music", nrow(dfp.wide))
                         ))
proportionsBulk <- rbind(proportionsBulk,
                         data.frame(
                           glial_known = 1-dfp.wide$neuron.true,
                           neuron_known = dfp.wide$neuron.true,
                           glial__predicted = 1-dfp.wide$neuron.music.noscale,
                           neuron_predicted = dfp.wide$neuron.music.noscale,
                           sample.id = rownames(dfp.wide),
                           condition.id = rep("no", nrow(dfp.wide)),
                           replicate.id = rep("music", nrow(dfp.wide))
                         ))
# bisque
proportionsBulk <- rbind(proportionsBulk,
                         data.frame(
                           glial_known = 1-dfp.wide$neuron.true,
                           neuron_known = dfp.wide$neuron.true,
                           glial__predicted = 1-dfp.wide$neuron.bisque.scale,
                           neuron_predicted = dfp.wide$neuron.bisque.scale,
                           sample.id = rownames(dfp.wide),
                           condition.id = rep("scale", nrow(dfp.wide)),
                           replicate.id = rep("bisque", nrow(dfp.wide))
                         ))
proportionsBulk <- rbind(proportionsBulk,
                         data.frame(
                           glial_known = 1-dfp.wide$neuron.true,
                           neuron_known = dfp.wide$neuron.true,
                           glial__predicted = 1-dfp.wide$neuron.bisque.noscale,
                           neuron_predicted = dfp.wide$neuron.bisque.noscale,
                           sample.id = rownames(dfp.wide),
                           condition.id = rep("no", nrow(dfp.wide)),
                           replicate.id = rep("bisque", nrow(dfp.wide))
                         ))

# Cohort 3 pseudobulk
load(file.path("./monaco-et-al-2019_pbmc-cohort1/env/01_pseudobulk/04_simulation.RData"))
proportionsCohort3 <- dfPseudobulk
# 
colnames(proportionsCohort3) <- c(
  "nonplasmablasts_predicted", "plasmablasts_predicted",
  "condition.id", "plasmablasts_known", "nonplasmablasts_known",
  "sample.id"
)
proportionsCohort3 <- proportionsCohort3[,c(1,2,4,5,6,3)]

# Cohort 1, proportions shuffle
load(file.path("./huuki-meyers-et-al-2023_dlpfc-cohort1/env/03_shuffle/02_figde_script.RData"))

proportionsShuffle <- dfp.tall
colnames(proportionsShuffle) <- c(
  "glial_known", "neuron_known", "glial_predicted", "neuron_predicted",
  "sample.id.shuffle.pred", "condition.id", "sample.id", "matched.id"
)
proportionsShuffle <- proportionsShuffle[,c(1,2,3,4,7,6)]
rm(dfp.tall)

# Run RMSE calculations
# K2 calculations
# cohorts 1 and 2,
listProportionsTablesK2 <- list(
  proportionsCohort1K2, proportionsCohort2K2, proportionsShuffle
)

# 2 cell types
experimentLabelsVectorK2 <- 
  c("huuki-meyers-et-al-2023_dlpfc-cohort1_k2", 
    "tran-et-al-2021_dlpfc-cohort2_k2", 
    "huuki-meyers-et-al-2023_dlpfc-cohort1_shuffle")
cellTypesVectorK2 <- c("glial", "neuron")
rmseSuppTableK2Twotypes <- do.call(rbind, 
                                   lapply(seq(length(listProportionsTablesK2)),
                                          function(index){
                                            message(index)
                                            proportionsData <- listProportionsTablesK2[[index]]
                                            if(!"replicate.id" %in% colnames(proportionsData)){
                                              proportionsData$replicate.id <- "rep1"
                                            }
                                            dataLabel <- experimentLabelsVectorK2[[index]]
                                            returnList <- applyRmse(
                                              proportionsData, 
                                              cellTypesVector=cellTypesVectorK2
                                            )
                                            returnTable <- returnList$returnRmseTall
                                            returnTable$experimentLabel <- experimentLabelsVectorK2[index]
                                            return(returnTable)
                                          })) |> 
  as.data.frame()

# 1 cell types
experimentLabelsVectorK2 <- 
  c("huuki-meyers-et-al-2023_dlpfc-cohort1_neuron_k2", 
    "tran-et-al-2021_dlpfc-cohort2_neuron_k2", 
    "huuki-meyers-et-al-2023_dlpfc-cohort1_neuron_shuffle")
cellTypesVectorK2 <- c("neuron")
rmseSuppTableK2Neuron <- do.call(rbind, 
                                 lapply(seq(length(listProportionsTablesK2)),
                                        function(index){
                                          message(index)
                                          proportionsData <- listProportionsTablesK2[[index]]
                                          if(!"replicate.id" %in% colnames(proportionsData)){
                                            proportionsData$replicate.id <- "rep1"
                                          }
                                          dataLabel <- experimentLabelsVectorK2[[index]]
                                          returnList <- applyRmse(
                                            proportionsData, 
                                            cellTypesVector=c("neuron")
                                          )
                                          returnTable <- returnList$returnRmseTall
                                          returnTable$experimentLabel <- experimentLabelsVectorK2[index]
                                          return(returnTable)
                                        })) |> 
  as.data.frame()
experimentLabelsVectorK2 <- 
  c("huuki-meyers-et-al-2023_dlpfc-cohort1_glial_k2", 
    "tran-et-al-2021_dlpfc-cohort2_glial_k2", 
    "huuki-meyers-et-al-2023_dlpfc-cohort1_glial_shuffle")
cellTypesVectorK2 <- c("glial")
rmseSuppTableK2Glial <- do.call(rbind, 
                                lapply(seq(length(listProportionsTablesK2)),
                                       function(index){
                                         message(index)
                                         proportionsData <- listProportionsTablesK2[[index]]
                                         if(!"replicate.id" %in% colnames(proportionsData)){
                                           proportionsData$replicate.id <- "rep1"
                                         }
                                         dataLabel <- experimentLabelsVectorK2[[index]]
                                         returnList <- applyRmse(
                                           proportionsData, 
                                           cellTypesVector=c("glial")
                                         )
                                         returnTable <- returnList$returnRmseTall
                                         returnTable$experimentLabel <- experimentLabelsVectorK2[index]
                                         return(returnTable)
                                       })) |> 
  as.data.frame()

# cohort 3, 
# 2 cell types
cellTypesVectorK2 <- c("nonplasmablasts", "plasmablasts")
experimentLabelsVectorK2 <- c("monaco-et-al-2019_pbmc-cohort1_k2")
listProportionsTablesK2 <- list(
  proportionsCohort3
)
rmseSuppTableK2Cohort3Twotypes <- do.call(rbind, 
                                          lapply(seq(length(listProportionsTablesK2)),
                                                 function(index){
                                                   message(index)
                                                   proportionsData <- listProportionsTablesK2[[index]]
                                                   if(!"replicate.id" %in% colnames(proportionsData)){
                                                     proportionsData$replicate.id <- "rep1"
                                                   }
                                                   dataLabel <- experimentLabelsVectorK2[[index]]
                                                   returnList <- applyRmse(
                                                     proportionsData, 
                                                     cellTypesVector=cellTypesVectorK2
                                                   )
                                                   returnTable <- returnList$returnRmseTall
                                                   returnTable$experimentLabel <- experimentLabelsVectorK2[index]
                                                   return(returnTable)
                                                 })) |> 
  as.data.frame()

# 1 cell types
cellTypesVectorK2 <- c("nonplasmablasts")
experimentLabelsVectorK2 <- c("monaco-et-al-2019_pbmc-cohort1_nonplasmablasts_k2")
listProportionsTablesK2 <- list(
  proportionsCohort3
)
rmseSuppTableK2Cohort3Nonplasmablasts <- do.call(rbind, 
                                                 lapply(seq(length(listProportionsTablesK2)),
                                                        function(index){
                                                          message(index)
                                                          proportionsData <- listProportionsTablesK2[[index]]
                                                          if(!"replicate.id" %in% colnames(proportionsData)){
                                                            proportionsData$replicate.id <- "rep1"
                                                          }
                                                          dataLabel <- experimentLabelsVectorK2[[index]]
                                                          returnList <- applyRmse(
                                                            proportionsData, 
                                                            cellTypesVector=cellTypesVectorK2
                                                          )
                                                          returnTable <- returnList$returnRmseTall
                                                          returnTable$experimentLabel <- experimentLabelsVectorK2[index]
                                                          return(returnTable)
                                                        })) |> 
  as.data.frame()

cellTypesVectorK2 <- c("plasmablasts")
experimentLabelsVectorK2 <- c("monaco-et-al-2019_pbmc-cohort1_plasmablasts_k2")
listProportionsTablesK2 <- list(
  proportionsCohort3
)
rmseSuppTableK2Cohort3Plasmablasts <- do.call(rbind, 
                                              lapply(seq(length(listProportionsTablesK2)),
                                                     function(index){
                                                       message(index)
                                                       proportionsData <- listProportionsTablesK2[[index]]
                                                       if(!"replicate.id" %in% colnames(proportionsData)){
                                                         proportionsData$replicate.id <- "rep1"
                                                       }
                                                       dataLabel <- experimentLabelsVectorK2[[index]]
                                                       returnList <- applyRmse(
                                                         proportionsData, 
                                                         cellTypesVector=cellTypesVectorK2
                                                       )
                                                       returnTable <- returnList$returnRmseTall
                                                       returnTable$experimentLabel <- experimentLabelsVectorK2[index]
                                                       return(returnTable)
                                                     })) |> 
  as.data.frame()

# try k2 from bulk
cellTypesVectorK2 <- c("glial", "neuron")
experimentLabelsVectorK2 <- c(
  "huuki-meyers-et-al-2023_dlpfc-cohort1_k2_bulk_nnls", 
  "huuki-meyers-et-al-2023_dlpfc-cohort1_k2_bulk_music", 
  "huuki-meyers-et-al-2023_dlpfc-cohort1_k2_bulk_bisque")
listProportionsTablesK2 <- list(
  proportionsBulk[proportionsBulk$replicate.id=="nnls",],
  proportionsBulk[proportionsBulk$replicate.id=="music",],
  proportionsBulk[proportionsBulk$replicate.id=="bisque",]
)
rmseSuppTableK2Bulk <- do.call(rbind, 
                               lapply(seq(length(listProportionsTablesK2)),
                                      function(index){
                                        message(index)
                                        proportionsData <- listProportionsTablesK2[[index]]
                                        if(!"replicate.id" %in% colnames(proportionsData)){
                                          proportionsData$replicate.id <- "rep1"
                                        }
                                        dataLabel <- experimentLabelsVectorK2[[index]]
                                        returnList <- applyRmse(
                                          proportionsData, 
                                          cellTypesVector=cellTypesVectorK2
                                        )
                                        returnTable <- returnList$returnRmseTall
                                        returnTable$experimentLabel <- experimentLabelsVectorK2[index]
                                        return(returnTable)
                                      })) |> 
  as.data.frame()

# K3 calculations
# cohort 1 and cohort2

# 3 cell types
cellTypesVectorK3 <- c("glial", "Inhib", "Excit")
experimentLabelsVectorK3 <- 
  c("huuki-meyers-et-al-2023_k3", "tran-et-al-2021_dlpfc-cohort2_k3")
listProportionsTablesK3 <- list(
  proportionsCohort1K3, proportionsCohort2K3
)
rmseSuppTableK3Threetypes <- do.call(rbind, 
                                     lapply(seq(length(listProportionsTablesK3)),
                                            function(index){
                                              message(index)
                                              proportionsData <- listProportionsTablesK3[[index]]
                                              if(!"replicate.id" %in% colnames(proportionsData)){
                                                proportionsData$replicate.id <- "rep1"
                                              }
                                              dataLabel <- experimentLabelsVectorK3[[index]]
                                              returnList <- applyRmse(
                                                proportionsData, 
                                                cellTypesVector=cellTypesVectorK3
                                              )
                                              returnTable <- returnList$returnRmseTall
                                              returnTable$experimentLabel <- experimentLabelsVectorK3[index]
                                              return(returnTable)
                                            })) |> 
  as.data.frame()

# 1 cell types
cellTypesVectorK3 <- c("glial")
experimentLabelsVectorK3 <- 
  c("huuki-meyers-et-al-2023_glial_k3", "tran-et-al-2021_dlpfc-cohort2_glial_k3")
rmseSuppTableK3Glial <- do.call(rbind, 
                                lapply(seq(length(listProportionsTablesK3)),
                                       function(index){
                                         message(index)
                                         proportionsData <- listProportionsTablesK3[[index]]
                                         if(!"replicate.id" %in% colnames(proportionsData)){
                                           proportionsData$replicate.id <- "rep1"
                                         }
                                         dataLabel <- experimentLabelsVectorK3[[index]]
                                         returnList <- applyRmse(
                                           proportionsData, 
                                           cellTypesVector=cellTypesVectorK3
                                         )
                                         returnTable <- returnList$returnRmseTall
                                         returnTable$experimentLabel <- experimentLabelsVectorK3[index]
                                         return(returnTable)
                                       })) |> 
  as.data.frame()

cellTypesVectorK3 <- c("Inhib")
experimentLabelsVectorK3 <- 
  c("huuki-meyers-et-al-2023_Inhib_k3", "tran-et-al-2021_dlpfc-cohort2_Inhib_k3")
rmseSuppTableK3Inhib <- do.call(rbind, 
                                lapply(seq(length(listProportionsTablesK3)),
                                       function(index){
                                         message(index)
                                         proportionsData <- listProportionsTablesK3[[index]]
                                         if(!"replicate.id" %in% colnames(proportionsData)){
                                           proportionsData$replicate.id <- "rep1"
                                         }
                                         dataLabel <- experimentLabelsVectorK3[[index]]
                                         returnList <- applyRmse(
                                           proportionsData, 
                                           cellTypesVector=cellTypesVectorK3
                                         )
                                         returnTable <- returnList$returnRmseTall
                                         returnTable$experimentLabel <- experimentLabelsVectorK3[index]
                                         return(returnTable)
                                       })) |> 
  as.data.frame()

cellTypesVectorK3 <- c("Excit")
experimentLabelsVectorK3 <- 
  c("huuki-meyers-et-al-2023_Excit_k3", "tran-et-al-2021_dlpfc-cohort2_Excit_k3")
rmseSuppTableK3Excit <- do.call(rbind, 
                                lapply(seq(length(listProportionsTablesK3)),
                                       function(index){
                                         message(index)
                                         proportionsData <- listProportionsTablesK3[[index]]
                                         if(!"replicate.id" %in% colnames(proportionsData)){
                                           proportionsData$replicate.id <- "rep1"
                                         }
                                         dataLabel <- experimentLabelsVectorK3[[index]]
                                         returnList <- applyRmse(
                                           proportionsData, 
                                           cellTypesVector=cellTypesVectorK3
                                         )
                                         returnTable <- returnList$returnRmseTall
                                         returnTable$experimentLabel <- experimentLabelsVectorK3[index]
                                         return(returnTable)
                                       })) |> 
  as.data.frame()

# Get aggregated table
# Format values across tables on aggregation.

roundDigits <- 3

listRmseTables <- list(
  rmseSuppTableK2Twotypes,
  rmseSuppTableK2Neuron,
  rmseSuppTableK2Glial,
  rmseSuppTableK3Threetypes,
  rmseSuppTableK3Excit,
  rmseSuppTableK3Inhib,
  rmseSuppTableK3Glial,
  rmseSuppTableK2Bulk,
  rmseSuppTableK2Cohort3Twotypes,
  rmseSuppTableK2Cohort3Nonplasmablasts,
  rmseSuppTableK2Cohort3Plasmablasts
)

# format rmse values
listRmseTables <- lapply(
  listRmseTables, function(tableIteration){
    tableIteration$rmseResult <- format(
      tableIteration$rmseResult, digits = roundDigits,
      scientific = TRUE
    )
    return(tableIteration)
  }
)

#----------------
# get final table
#----------------
rmseSuppTableOut <- do.call(rbind, 
                            lapply(listRmseTables, function(item){item})) |> as.data.frame()

removeSampleRmse <- !grepl("sample.id", rmseSuppTableOut$condition)
table(removeSampleRmse)
rmseSuppTableOut <- rmseSuppTableOut[removeSampleRmse,]
head(rmseSuppTableOut)

rmseSuppTableOut <- rmseSuppTableOut[!grepl('rep1;.*', rmseSuppTableOut$condition),]
rmseSuppTableOut$condition <- gsub("replicate.id:rep1", "all", rmseSuppTableOut$condition)
rmseSuppTableOut$condition <- gsub("replicate.id:", "", rmseSuppTableOut$condition)
rmseSuppTableOut$condition <- gsub("condition.id:", "", rmseSuppTableOut$condition)
rmseSuppTableOut

rmseSuppTableOut$condition <- gsub("^(scaled|scale)$", "withscale", rmseSuppTableOut$condition)
rmseSuppTableOut$condition <- gsub(";(scaled|scale)$", ";withscale", rmseSuppTableOut$condition)
rmseSuppTableOut$condition <- gsub("^(no|unscaled)$", "noscale", rmseSuppTableOut$condition)
rmseSuppTableOut$condition <- gsub(";(no|unscaled)$", ";noscale", rmseSuppTableOut$condition)
rmseSuppTableOut

# Save final RMSE table
rmseSupplement <- rmseSuppTableOut
# append columns
# append cohort
rmseSupplement$cohort <- ifelse(
  grepl("huuki-meyers-et-al-2023", rmseSupplement$experimentLabel), "Huuki-Myers et al 2023",
  ifelse(grepl("tran-et-al-2021", rmseSupplement$experimentLabel), "Tran et al 2021",
         ifelse(grepl("monaco-et-al-2019", rmseSupplement$experimentLabel), "Monaco et al 2019", "NA")))
# append k
rmseSupplement$k <- ifelse(grepl("k2|shuffle", rmseSupplement$experimentLabel), "2",
                           ifelse(grepl("k3", rmseSupplement$experimentLabel), "3", "NA"))

# append cell type strings
rmseSupplement$cellTypes <- ""
# k2
filterK2Brain <- rmseSupplement$k == "2" & 
  grepl("huuki-meyers-et-al-2023|tran-et-al-2021", rmseSupplement$experimentLabel)
filterK2Brain2Types <- filterK2Brain & rmseSupplement$numCellTypes=="2"
rmseSupplement$cellTypes[filterK2Brain2Types] <- "neuron;glial"
filterK2BrainNeuron <- filterK2Brain & grepl("neuron", rmseSupplement$experimentLabel)
rmseSupplement$cellTypes[filterK2BrainNeuron] <- "neuron"
filterK2BrainGlial <- filterK2Brain & grepl("glial", rmseSupplement$experimentLabel)
rmseSupplement$cellTypes[filterK2BrainGlial] <- "glial"
# k3
filterK3Brain <- rmseSupplement$k == "3" & 
  grepl("huuki-meyers-et-al-2023|tran-et-al-2021", rmseSupplement$experimentLabel)
filterK3Brain3Types <- filterK3Brain & rmseSupplement$numCellTypes==3
rmseSupplement$cellTypes[filterK3Brain3Types] <- "Excit;Inhib;glial"
filterK3BrainExcit <- filterK3Brain & grepl("Excit", rmseSupplement$experimentLabel)
rmseSupplement$cellTypes[filterK3BrainExcit] <- "Excit"
filterK3BrainInhib <- filterK3Brain & grepl("Inhib", rmseSupplement$experimentLabel)
rmseSupplement$cellTypes[filterK3BrainInhib] <- "Inhib"
filterK3BrainGlial <- filterK3Brain & grepl("glial", rmseSupplement$experimentLabel)
rmseSupplement$cellTypes[filterK3BrainGlial] <- "glial"
# k2 pbmc
filterK2Blood <- grepl("monaco-et-al-2019", rmseSupplement$experimentLabel)
filterK2Blood2Types <- filterK2Blood & rmseSupplement$numCellTypes==2
rmseSupplement$cellTypes[filterK2Blood2Types] <- "plasmablasts;nonplasmablasts"
filterK2BloodPlasmablasts <- filterK2Blood & grepl("plasmablasts", rmseSupplement$experimentLabel)
rmseSupplement$cellTypes[filterK2BloodPlasmablasts] <- "plasmablasts"
filterK2BloodNonplasmablasts <- filterK2Blood & grepl("nonplasmablasts", rmseSupplement$experimentLabel)
rmseSupplement$cellTypes[filterK2BloodNonplasmablasts] <- "nonplasmablasts"

filterK3Brain <- rmseSupplement$k == "3" & grepl("cohort1|cohort2", rmseSupplement$experimentLabel)
rmseSupplement$cellTypes[filterK3Brain] <- "Excit;Inhib;glial"

# append num cell types
rmseSupplement$numCellTypes <- gsub("k", "", rmseSupplement$k)

singletonCellTypesVector <- c("neuron", "glial", "plasmablasts", "nonplasmablasts", "Inhib", "Excit")
singletonCellTypesVector <- paste0("^", singletonCellTypesVector, "$")
singletonCellTypesVector <- paste(singletonCellTypesVector, collapse = "|")

conditionCellTypes <- grepl(singletonCellTypesVector, rmseSupplement$cellTypes)
rmseSupplement$numCellTypes[conditionCellTypes] <- 1

# remove conditions
dim(rmseSupplement)
filterCondition <- grepl("sample.id", rmseSupplement$condition)
filterCondition <- filterCondition & grepl("cohort3", rmseSupplement$experimentLabel)
rmseSupplement <- rmseSupplement[!filterCondition,]
dim(rmseSupplement)

rmseSuppTableOutRound <- rmseSupplement
rmseSuppTableOutRound$rmseResult <- 
  format(rmseSuppTableOutRound$rmseResult, digits = 3)

#-----
# save
#-----
#
#
save(
  rmseSupplement, 
  file = "./software_manuscript/outputs/02_aggregateRMSE/rmseSupplement.rda"
  )
