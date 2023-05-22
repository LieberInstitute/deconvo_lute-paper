#!/usr/bin/env R

# Author: Sean Maden
#
# Joint analysis of halo, bulk, and sn data, and potential confounds
#
#

libv <- c("dplyr", "ggplot2")
sapply(libv, library, character.only = T)

#---------------------------------
# prepare processed confounds data
#---------------------------------
# load confounds
halo.confounds <- get(load(
  "./deconvo_method-paper/outputs/22_factors-impacting-deconvo-errors/list-halo-confounds_dlpfc-train.rda"))
bulk.confounds <- get(load(
  "./deconvo_method-paper/outputs/22_factors-impacting-deconvo-errors/list-bulk-confounds_dlpfc-train.rda"))
sn.confounds <- get(load(
  "./deconvo_method-paper/outputs/22_factors-impacting-deconvo-errors/df-confounds-sn_dlpfc-train.rda"))

# load marker expression
#bulk.marker.expr <- get(load("list-donormarker-bulkexpr_dlpfc-train.rda"))
#sn.marker.expr <- get(load("list-donormarker-snexpr_dlpfc-train.rda"))

#---------------------
# harmonize assay data
#---------------------
library(dplyr)

sample.id.varname <- "sample.id"
df.bulk <- cbind(bulk.confounds$total.counts, bulk.confounds$expressed.genes)

df.bulk[,sample.id.varname] <- paste0(
  unlist(strsplit(df.bulk$sample.id, "_"))[2], "_",
  tolower(unlist(strsplit(df.bulk$sample.id, "_"))[3]))

for(sample.id in rownames(df.sn.confound)){
  for(cn in colnames(df.sn.confound)){
    filter.df.sn <- df.sn.confound$sample.id==sample.id
    df.bulk[df.bulk$sample.id.sn==sample.id, cn] <- df.sn.confound[filter.df.sn, cn]
  }
}

df.halo <- do.call(cbind, halo.confounds)
df.halo <- df.halo[,c(1,2,4,6)]
colnames(df.halo) <- c("sample.id.halo", "total.cells", "total.neuron", "total.glial")
df.halo$sample.id <- gsub("_.*", "", df.halo[,1])
df.halo$sample.id <- paste0(gsub("M$|P$|A$", "", df.halo$sample.id), "_",
                            ifelse(grepl("M", df.halo$sample.id), "mid",
                                   ifelse(grepl("A", df.halo$sample.id), "ant", "post")))
df.halo.agg <- aggregate(x = list(df.halo$total.cells, df.halo$total.neuron, 
                                  df.halo$total.glial), 
                         by = list(df.halo$sample.id), FUN = "mean")
colnames(df.halo.agg) <- c("sample.id", "halo.total.cells", 
                           "halo.total.neuron", "halo.total.glial")

for(sample.id in df.halo.agg$sample.id){
  for(cn in colnames(df.halo.agg)){
    filter.halo <- df.halo.agg$sample.id==sample.id
    df.bulk[df.bulk$sample.id==sample.id,cn] <- df.halo.agg[filter.halo,cn]
  }
}

# format integrated data
df.integrate <- df.bulk
df.integrate$cell.compartment <- gsub(".*_", "", df.integrate$sample.id.1)
df.integrate$library.prep <- gsub("_.*", "", df.integrate$sample.id.1)

# save
df.integrate.path <- paste0("./deconvo_method-paper/outputs/",
                            "22_factors-impacting-deconvo-errors/",
                            "df_integrate.rda")
save(df.integrate, file = df.integrate.path)

#----------------------
# cross-assay summaries
#----------------------
df.integrate.path <- paste0("./deconvo_method-paper/outputs/",
                            "22_factors-impacting-deconvo-errors/",
                            "df_integrate.rda")
df.integrate <- get(load(df.integrate.path))

# plot summaries

ggplot(df.integrate, aes(x = sample.id, y = expressed.genes)) + 
  geom_violin(draw_quantiles = 0.5)

ggplot(df.integrate, aes(x = halo.total.glial, y = expressed.genes)) + 
  geom_point(alpha = 0.5) + geom_abline(slope = 1, intercept = 0)




