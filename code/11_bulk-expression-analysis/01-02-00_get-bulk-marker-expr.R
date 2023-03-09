#!/usr/bin/env R

# Author: Sean Maden
#
#

libv <- c("SingleCellExperiment", "SummarizedExperiment")
sapply(libv, library, character.only = T)

#----------
# load data
#----------
# load bulk data
rse.filename <- "rse_gene.Rdata"
load.path <- file.path("Human_DLPFC_Deconvolution", "processed-data",
                       "01_SPEAQeasy", "round2_v40_2022-07-06", "rse")
rse <- get(load(file.path(load.path, rse.filename)))

# get save directory path
save.path <- file.path("deconvo_method-paper", "outputs", 
                       "11_bulk-expression-analysis")

# load marker data
sce.markers.path <- file.path("deconvo_method-paper", "outputs", "09_manuscript", 
                              "list-scef_markers-k2-k3-k4_ro1-dlpfc.rda")
lsce <- get(load(sce.markers.path))
sce <- lsce[["k2"]]
rm(lsce)

# load helper functions
group.jitter.filename <- "group-jitter_helper.rda"
get.summary.filename <- "get-summary_helper.rda"
group_jitter <- get(load(file.path(save.path, group.jitter.filename)))
get_summary_list <- get(load(file.path(save.path, get.summary.filename)))

#------------------------
# manage helper functions
#------------------------

ggplot_from_dfp <- function(dfp, variable.name,
                            category.type = c("group", "color", "fill"),
                            plot.type = c("box", "jitter", "both"),
                            scale = c("normal", "log")){
  # make new plot object
  message("Getting new plot object...")
  
  # parse plot formatting
  if(category.type == "group"){
    new.plot <- ggplot(dfp, aes(x = group, y = value, group = category))
  } else if(category.type == "color"){
    new.plot <- ggplot(dfp, aes(x = group, y = value, color = category))
  } else if(category == "fill"){
    new.plot <- ggplot(dfp, aes(x = group, y = value, fill = category))
  } else if(category == "fill;color"){
    new.plot <- ggplot(dfp, aes(x = group, y = value, 
                                fill = category, color = category))
  } else{}
  
  # parse plot options
  if(plot.type == "box"){
    new.plot <- new.plot + geom_boxplot()
  } else if(plot.type == "jitter"){
    new.plot <- new.plot + geom_jitter()
  } else{
    new.plot <- new.plot + geom_jitter() + geom_boxplot()
  }
  
  if(scale == "log"){new.plot <- new.plot + scale_y_log10()}
  
  # parse theme and axes
  new.plot <- new.plot + theme_bw() + ylab(variable.name) + xlab(variable.name) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
  return(new.plot)
}

get_comparison_data <- function(expr.bg, expr.marker, plot.fname,
                                type.vector, cd, save.path){
  for(type in type.vector){
    message("Working on summary type ", type, "...")
    message("Working on background expression...")
    lbg <- get_summary_list(type = type, plot.fname = plot.fname, 
                            variable.vector = variable.vector, cd = cd, 
                            counts = expr.bg, save.path = save.path)
    message("Working on marker expression...")
    lmarker <- get_summary_list(type = type, plot.fname = plot.fname, 
                                variable.vector = variable.vector, cd = cd, 
                                counts = expr.marker, save.path = save.path)
    message("Getting new expression data...")
    ldfp <- lapply(seq(length(lbg)), function(i){
      dfp.bg <- lbg[[i]][[1]]; dfp.marker <- lmarker[[i]][[1]]
      dfp.bg$marker.type <- "background"; dfp.marker$marker.type <- "marker"
      return(rbind(dfp.bg, dfp.marker))
    })
    message("Finished with summary type ", type, ".")
    
  }
  message("Done with all summary  types.")
  message("Getting plots...")
  # get boxplots
  lbox <- lapply(seq(length(ldfp)), function(i){
    dfp <- ldfp[[i]]; dfp$category <- dfp$marker.type
    ggplot_from_dfp(dfp, type.vector[i], "color", "box", scale = "normal")
  })
  # get boxplots -- log scale yaxis
  lbox.logscale <- lapply(seq(length(ldfp)), function(i){
    dfp <- ldfp[[i]]; dfp$category <- dfp$marker.type
    ggplot_from_dfp(dfp, type.vector[i], "color", "box", scale = "log")
  })
  message("Getting return list...")
  lr <- list(dfp = ldfp, lgg = list(box = lbox, box.log = lbox.logscale))
  return(lr)
}

#--------------------------------
# get marker gene bulk expression
#--------------------------------
# subset rse
bulk.gene.names <- rowData(rse)$Symbol
marker.genes.vector <- rownames(sce)
overlapping.markers <- intersect(bulk.gene.names, marker.genes.vector)
message("Found ", length(overlapping.markers), " overlapping markers.")
filter <- which(rowData(rse)$Symbol %in% overlapping.markers)
rsef <- rse[filter,]
dim(rsef)

# save bulk marker expr
rsef.filename <- "rsef_k2-marker-expr_ro1-dlpfc.rda"
save.path <- file.path(save.path, rsef.filename)
save(rsef, file = save.path)

#--------------------------
# prepare comparison params
#--------------------------
# set up data summaries
# params
assay.name <- "counts"
batch.variable <- "batch.id"
condition.variable <- "expt_condition"
# rse data
cd <- colData(rse)
# get new cd variables
cd[,batch.variable] <- paste0(cd$BrNum,"_",cd$location)
cd[,condition.variable] <- paste0(cd$library_prep,"_",cd$library_type)
# vector of group variables to summarize
variable.vector <- c(batch.variable, "library_prep", "library_type", 
                     condition.variable)
# define summary types
type.vector <- c("total.counts", "zero.count", "mean", "variance")

#---------------
# compare counts
#---------------
counts.bg <- assays(rse)[["counts"]]
counts.marker <- assays(rsef)[["counts"]]

lcomp <- get_comparison_data(expr.bg = counts.bg, 
                             expr.marker = counts.marker, 
                             plot.fname = plot.fname,
                             type.vector = type.vector, 
                             cd = cd, save.path = save.path)


#--------------------------
# compare normalized counts
#--------------------------

#------------
# compare tpm
#------------

#-------------
# compare fpkm
#-------------

#--------------------------------
# qc at markers versus background
#--------------------------------








