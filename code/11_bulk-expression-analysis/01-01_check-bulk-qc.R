#!/usr/bin/env R

# Author: Sean Maden
#
#

libv <- c("SummarizedExperiment", "ggplot2")
sapply(libv, library, character.only = T)

#----------
# load data
#----------
# load bulk data
rse.filename <- "rse_gene.Rdata"
load.path <- file.path("Human_DLPFC_Deconvolution",
                       "processed-data",
                       "01_SPEAQeasy",
                       "round2_v40_2022-07-06",
                       "rse")
rse <- get(load(file.path(path, rse.filename)))

# get save directory path
save.path <- file.path("deconvo_method-paper", "outputs", 
                       "11_bulk-expression-analysis")

#---------------------------
# total expression summaries
#---------------------------
# params
assay.name <- "counts"
batch.variable <- "batch.id"
# rse data
cd <- colData(rse)
cd[,batch.variable] <- paste0(cd$BrNum,"_",cd$location)
counts <- assays(rse)[[assay.name]]

# total expression by sample
dfp <- data.frame(total.counts = colSums(counts))
dfp$sample <- cd[,batch.variable]
labels.vector <- unique(dfp$sample)
levels.vector <- unlist(sapply(labels.vector, function(li){median(dfp[dfp[,2]==li,1])}))
dfp$sample <- factor(dfp$sample, levels = labels.vector[order(levels.vector)])
ggplot(dfp, aes(x = sample, y = total.counts)) + theme_bw() +
  geom_jitter() + geom_boxplot(color = "cyan", alpha = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# total expression by library prep

# total expression by condition

# total expression by condition, library prep


#-----------------------------------
# missing/unexpressed gene summaries
#-----------------------------------

#------------------------
# mean-variance summaries
#------------------------
