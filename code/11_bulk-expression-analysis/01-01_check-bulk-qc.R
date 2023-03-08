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

#-----------------
# helper functions
#-----------------
group_jitter <- function(variable.name, cd, counts,
                         ggtitle.string = "Total counts by group",
                         type = "total.counts"){
  # get plot data
  group.vector <- unique(cd[,variable.name])
  if(type == "total.counts"){
    value.string <- "Total counts"
    dfp <- do.call(rbind, lapply(group.vector, function(gi){
      filter <- cd[,variable.name]==gi
      dfi <- data.frame(value = colSums(counts[,filter]))
      dfi$group <- gi; return(dfi)
    }))
  } else if(type == "number.zero"){
    value.string <- "Num. zero"
    dfp <- do.call(rbind, lapply(group.vector, function(gi){
      filter <- cd[,variable.name]==gi
      cf <- counts[,filter]
      
      num.zero <- apply(cf,2,function(ci){length(ci[ci==0])})
      
      dfi <- data.frame(value = colSums(counts[,filter]))
      dfi$group <- gi; return(dfi)
    }))
  } else{
    stop("Error, unrecognized type argument.")
  }
  
  
  # order group as factor
  labels.vector <- unique(dfp$group)
  levels.vector <- unlist(sapply(labels.vector, function(li){
    median(dfp[dfp[,2]==li,1])}))
  dfp$group <- factor(dfp$group, 
                      levels = labels.vector[order(levels.vector)])
  
  # make new plot object
  new.plot <- ggplot(dfp, aes(x = group, y = total.counts)) + theme_bw() +
    geom_jitter() + geom_boxplot(color = "cyan", alpha = 0) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(ggtitle.string) + xlab(value.string)
  return(list(dfp = dfp, new.plot = new.plot))
}


#---------------------------
# total expression summaries
#---------------------------
# params
assay.name <- "counts"
batch.variable <- "batch.id"
condition.variable <- "expt_condition"
# rse data
cd <- colData(rse)
counts <- assays(rse)[[assay.name]]
# get new cd variables
cd[,batch.variable] <- paste0(cd$BrNum,"_",cd$location)
cd[,condition.variable] <- paste0(cd$library_prep,"_",cd$library_type)

# total expression by sample
ljitter1 <- group_jitter(batch.variable, cd, counts, "Batch ID")

# total expression by library prep
ljitter2 <- group_jitter("library_prep", cd, counts, "Batch ID")

# total expression by condition
ljitter3 <- group_jitter("library_type", cd, counts, "Batch ID")

# total expression by condition, library prep
ljitter4 <- group_jitter(condition.variable, cd, counts, "Batch ID")

#-----------------------------------
# missing/unexpressed gene summaries
#-----------------------------------

#------------------------
# mean-variance summaries
#------------------------
