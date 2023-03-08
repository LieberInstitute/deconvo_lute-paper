#!/usr/bin/env R

# Author: Sean Maden
#
#

libv <- c("SummarizedExperiment", "ggplot2", "glmGamPoi", "gridExtra")
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
                         type = c("total.counts", "zero.count",
                                  "mean", "variance", "dispersion")){
  # define summary groups
  group.vector <- unique(cd[,variable.name])
  # get plot data
  message("Calculating summary data...")
  if(type == "total.counts"){
    value.string <- "Total counts"
    dfp <- do.call(rbind, lapply(group.vector, function(gi){
      message("Getting value ", type, " for group ", gi, "...")
      filter <- cd[,variable.name]==gi
      dfi <- data.frame(value = colSums(counts[,filter]))
      dfi$group <- gi; return(dfi)
    }))
  } else if(type == "zero.count"){
    value.string <- "Zero count"
    dfp <- do.call(rbind, lapply(group.vector, function(gi){
      message("Getting value ", type, " for group ", gi, "...")
      filter <- cd[,variable.name]==gi; cf <- counts[,filter]
      dfi <- data.frame(value = apply(cf,2,function(ci){
        length(which(ci==0))}))
      dfi$group <- gi; return(dfi)
    }))
  } else if(type == "variance"){
    value.string <- "Variance"
    dfp <- do.call(rbind, lapply(group.vector, function(gi){
      message("Getting value ", type, " for group ", gi, "...")
      filter <- cd[,variable.name]==gi
      dfi <- data.frame(value = colVars(counts[,filter]))
      dfi$group <- gi; return(dfi)
    }))
  } else if(type == "mean"){
    value.string <- "Mean"
    dfp <- do.call(rbind, lapply(group.vector, function(gi){
      message("Getting value ", type, " for group ", gi, "...")
      filter <- cd[,variable.name]==gi
      dfi <- data.frame(value = colMeans(counts[,filter]))
      dfi$group <- gi; return(dfi)
    }))
  } else if(type == "dispersion"){
    value.string <- "Dispersion"
    set.seed(0)
    num.gene <- 500
    cf <- counts[sample(seq(nrow(counts)), num.gene),]
    dfp <- do.call(rbind, lapply(group.vector, function(gi){
      message("Getting value ", type, " for group ", gi, "...")
      filter <- cd[,variable.name]==gi; cff <- cf[,filter]
      dispersion.vector <- apply(cff, 1, function(ri){
        glmGamPoi::glm_gp(ri)$overdispersions})
      dispersion.vector <- as.numeric(dispersion.vector)
      dfi <- data.frame(value = dispersion.vector)
      dfi$group <- gi; return(dfi)
    }))
  } else{
    stop("Error, unrecognized type argument.")
  }
  
  message("Formatting plot data...")
  # order group as factor
  labels.vector <- unique(dfp$group)
  levels.vector <- unlist(sapply(labels.vector, function(li){
    median(dfp[dfp[,2]==li,1])}))
  dfp$group <- factor(dfp$group, 
                      levels = labels.vector[order(levels.vector)])
  
  # make new plot object
  message("Getting new plot object...")
  new.plot <- ggplot(dfp, aes(x = group, y = value)) + theme_bw() +
    geom_jitter() + geom_boxplot(color = "cyan", alpha = 0) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab(value.string) + xlab(variable.name)
  
  return(list(dfp = dfp, new.plot = new.plot))
}

get_summary_list <- function(type, plot.fname, variable.vector, cd, counts, save.path){
  # get plot data
  
  ljitter <- lapply(variable.vector, function(variable){
    group_jitter(variable, cd, counts, type)
  })
  
  # save composite plot
  plot.fname <-  paste0("ggjitter-composite_", gsub("\\.", "-", type), 
                        "_ro1-dlpfc.jpg")
  plot.path <- file.path(save.path, plot.fname)
  message("Saving new plot: ", plot.path, "...")
  jpeg(plot.path, width = 10, height = 8, units = "in", res = 400)
  grid.arrange(ljitter[[1]][[2]], ljitter[[2]][[2]], ljitter[[3]][[2]], ljitter[[4]][[2]],
               nrow = 2, bottom = "group", left = "value")
  dev.off()
  return(ljitter)
}

#----------------------
# set up data summaries
#----------------------
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
# vector of group variables to summarize
variable.vector <- c(batch.variable, "library_prep", "library_type", 
                     condition.variable)

#----------------------
# save new summary data
#----------------------
# define summary types
type.vector <- c("total.counts", "dispersion", "zero.count", "mean", "variance")
for(type in type.vector){
  message("Working on summary type ", type, "...")
  get_summary_list(type = type, plot.fname = plot.fname, 
                   variable.vector = variable.vector, cd = cd, 
                   counts = counts, save.path = save.path)
  message("Finished with summary type ", type, ".")
}
message("Done with all summary  types.")
