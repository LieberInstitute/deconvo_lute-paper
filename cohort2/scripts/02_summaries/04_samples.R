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

sampleSummaries <- data.frame(
  samples = length(unique(cd$donor)),
  precent_female = sum(cd$sex == "F")/length(cd$sex),
  number_regions = length(unique(cd$region)),
  percent_posterior = sum(cd$region == "posterior")/length(cd$region),
  percent_middle = sum(cd$region == "middle")/length(cd$region),
  percent_anterior = sum(cd$region == "anterior")/length(cd$region),
  mean_nuclei_per_sample = mean(table(cd$donor)) |> round(0),
  median_nuclei_per_sample = median(table(cd$donor)) |> round(0),
  sd_nuclei_per_sample = sd(table(cd$donor)) |> round(0),
  total_nuclei = sum(table(cd$donor)) |> round(0)
)
sampleSummaries$percent_posterior <- 100

#-----
# save
#-----
write.csv(
  sampleSummaries, 
  file = "./outputs/02_summaries/cohort2_sample_summaries.csv", 
  row.names = FALSE)
