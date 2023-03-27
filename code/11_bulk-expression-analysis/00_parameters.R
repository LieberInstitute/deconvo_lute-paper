#!/usr/bin/env R

# Author: Sean Maden
#
# Parameters and depenencies for DLPFC RO1 bulk RNA-seq analyses.
#


# libraries
libv <- c("here", "SummarizedExperiment", "SingleCellExperiment", "dplyr", "data.table", "DESeq2", "ggplot2", "ggcorrplot", "glmGamPoi", "gridExtra")
sapply(libv, library, character.only = T)
# save path, outputs
save.path <- here("deconvo_method-paper", "outputs", "11_bulk-expression-analysis")

#-----------------
# helper functions
#-----------------
# summary/qc plots
group_jitter <- function(variable, cd, expression,
                         type = c("total.counts", "zero.count",
                                  "mean", "variance", "dispersion")){
  # define summary groups
  group.vector <- unique(cd[,variable])
  # get plot data
  message("Calculating summary data...")
  if(type == "total.counts"){
    value.string <- "Total counts"
    dfp <- do.call(rbind, lapply(group.vector, function(gi){
      message("Getting value ", type, " for group ", gi, "...")
      filter <- cd[,variable]==gi
      dfi <- data.frame(value = colSums(expression[,filter]))
      dfi$group <- gi; return(dfi)
    }))
  } else if(type == "zero.count"){
    value.string <- "Zero count"
    dfp <- do.call(rbind, lapply(group.vector, function(gi){
      message("Getting value ", type, " for group ", gi, "...")
      filter <- cd[,variable]==gi; cf <- expression[,filter]
      dfi <- data.frame(value = apply(cf,2,function(ci){
        length(which(ci==0))}))
      dfi$group <- gi; return(dfi)
    }))
  } else if(type == "variance"){
    value.string <- "Variance"
    dfp <- do.call(rbind, lapply(group.vector, function(gi){
      message("Getting value ", type, " for group ", gi, "...")
      filter <- cd[,variable]==gi
      dfi <- data.frame(value = colVars(expression[,filter]))
      dfi$group <- gi; return(dfi)
    }))
  } else if(type == "mean"){
    value.string <- "Mean"
    dfp <- do.call(rbind, lapply(group.vector, function(gi){
      message("Getting value ", type, " for group ", gi, "...")
      filter <- cd[,variable]==gi
      dfi <- data.frame(value = colMeans(expression[,filter]))
      dfi$group <- gi; return(dfi)
    }))
  } else if(type == "dispersion"){
    value.string <- "Dispersion"
    set.seed(0)
    num.gene <- 500
    # filter on protein-coding
    cf <- expression[sample(seq(nrow(expression)), num.gene),]
    dfp <- do.call(rbind, lapply(group.vector, function(gi){
      message("Getting value ", type, " for group ", gi, "...")
      filter <- cd[,variable]==gi; cff <- cf[,filter]
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
    ylab(value.string) + xlab(variable)
  
  return(list(dfp = dfp, new.plot = new.plot))
}

get_summary_list <- function(type, plot.filename, variable.vector, cd, expression, save.path){
  # get plot data
  
  ljitter <- lapply(variable.vector, function(variable){
    group_jitter(variable, cd, expression, type)
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

# bulk k2 marker expression
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
    new.plot <- ggplot(dfp, aes(x = group, y = value, fill = category, color = category))
  } else{}
  
  # parse plot options
  if(plot.type == "box"){
    new.plot <- new.plot + geom_boxplot()
  } else if(plot.type == "jitter"){
    new.plot <- new.plot + geom_jitter()
  } else{
    new.plot <- new.plot + geom_jitter() + geom_boxplot()
  }
  
  # parse theme and axes
  if(scale == "log"){
    new.plot <- new.plot + scale_y_log10() + theme_bw() + 
      ylab(paste0(variable.name, " (log-scale)")) + xlab(variable.name) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  } else{
    new.plot <- new.plot + theme_bw() + ylab(variable.name) + xlab(variable.name) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  return(new.plot)
}

get_comparison_data <- function(expr.bg, expr.marker, plot.filename,
                                type.vector, cd, save.path){
  for(type in type.vector){
    message("Working on summary type ", type, "...")
    message("Working on background expression...")
    lbg <- get_summary_list(type = type, plot.filename = plot.filename, 
                            variable.vector = variable.vector, cd = cd, 
                            expression = expr.bg, save.path = save.path)
    message("Working on marker expression...")
    lmarker <- get_summary_list(type = type, plot.filename = plot.filename, 
                                variable.vector = variable.vector, cd = cd, 
                                expression = expr.marker, save.path = save.path)
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

# bulk-pseudobulk marker expression
get_pseudobulk <- function(sce, sce.assay.name, S){
  sce.expression <- assays(sce)[[sce.assay.name]]
  sce.expression <- as.matrix(sce.expression)
  ZP <- do.call(cbind, lapply(K, function(ki){
    ki.filter <- colData(sce)[,sce.celltype.variable]==ki
    rowMeans(sce.expression[,ki.filter])
  }))
  colnames(ZP) <- K
  S <- S[order(match(names(S), colnames(ZP)))]
  if(!identical(colnames(ZP), names(S))){stop("Error, couldn't match cell type names in ZP, S.")}
  ZPS <- sweep(ZP, MARGIN = 2, STATS = S, FUN = "*")
  pseudobulk <- matrix(rowMeans(ZPS), ncol = 1)
  rownames(pseudobulk) <- rownames(sce)
  return(pseudobulk)
}

#------------
# params (01)
#------------
# gene type variables
gene.types.protein <- c("protein_coding")
gene.types.nonpolya <- c("lncRNA", "Mt_rRNA", "rRNA", "Mt_tRNA")
gene.types.include <- c(gene.types.protein, gene.types.nonpolya)

# paths
rse.bulk.filename <- "rse_gene.Rdata"
rse.bulk.filepath <- here("Human_DLPFC_Deconvolution", "processed-data", "01_SPEAQeasy", "round2_v40_2022-07-06", "rse", rse.bulk.filename)
rse.bulk.filename.new <- "rse-bulk-resave.rda"
rse.bulk.path.new <- here(save.path, rse.bulk.filename.new)
rse.k2markers.filename <- "rse_k2-marker-expression_ro1-dlpfc.rda"
rse.k2markers.filepath <- here(save.path, rse.k2markers.filename)
rse.gene.filter.filename <- "rse-gene-filter.rda"
rse.gene.filter.filepath <- here(save.path, rse.gene.filter.filename)
sce.markers.filename <- "list-scef_markers-k2-k3-k4_ro1-dlpfc.rda"
sce.markers.list.path <- here("deconvo_method-paper", "outputs", "09_manuscript", sce.markers.filename)

# variable names
experiment.condition1 <- "library_prep"
experiment.condition2 <- "library_type"
condition.variable <- "expt_condition"
donor.variable <- "BrNum"
location.variable <- "location"
batch.variable <- "batch.id"


#----------------------------------
# params for marker expression (05)
#----------------------------------
assay.name <- "counts"
assays <- c("counts")
batch.variable <- "batch.id"
condition.variable <- "expt_condition"
# vector of group variables to summarize
variable.vector <- c(batch.variable, "library_prep", "library_type", condition.variable)
# define summary types
type.vector <- c("total.counts", "zero.count", "mean", "variance")
correlation.method.markers <- "spearman"

#-----------------------------
# params for qc summaries (03)
#-----------------------------
# params
donor.id <- "BrNum"
brain.location.id <- "location"
assay.name.qc <- "counts"
batch.variable.qc <- "batch.id"
condition.variable.qc <- "expt_condition"
# vector of group variables to summarize
variable.vector <- c(batch.variable, "library_prep", "library_type", condition.variable)
type.vector.qc <- c("total.counts", "dispersion", "zero.count", "mean", "variance")

#----------------------------------
# params for marker expression (04)
#----------------------------------
condition.variable.qc <- "expt_condition"
sce.celltype.variable <- "k2"
sce.assay.name <- "counts"
K <- c("neuron", "glial")
S <- c("neuron" = 10,"glial" = 3)
correlation.heatmap.jpg.file.name <- "ggcorrheatmap_k2markers-rsegroup-pseudobulk.jpg"
correlation.heatmap.jpg.path <- here(save.path, correlation.heatmap.jpg.file.name)
pseudobulk.file.name <- "sce-pseudobulk.rda"
pseudobulk.path <- here(save.path, pseudobulk.file.name)

#-------------------------------------------------
# params for marker expression vs. background (05)
#-------------------------------------------------

#------------------------------------
# compare-k2markers-bulk-vs-halo (08)
#------------------------------------
# load halo outputs
halo.outputs.filename <- "halo_all.Rdata"
halo.outputs.path <- here("Human_DLPFC_Deconvolution", "processed-data", "03_HALO", halo.outputs.filename)
# halo table key

#-----------------------------
# compare mrna subspecies (09)
#-----------------------------
# histone rna filter
histone.patterns <- c("H1-", "H2A", "H2B", "H3C", "H3-", "H3P", "H3Y", "H4C")
histone.genes.pattern <- paste0("^", histone.patterns, ".*", collapse = "|")
# mitochondrial rna filter
mito.pattern <- c("MT-")
mito.genes.pattern <- paste0("^", mito.pattern, collapse = "|")
# new plot names
bulk.mean.rna.types.jpg.name <- "ggplot-jitter-boxplot_mean-expression_bulk-rna-types.jpg"
bulk.variance.rna.types.jpg.name <- "ggplot-jitter-boxplot_variance-expression_bulk-rna-types.jpg"
bulk.total.rna.types.jpg.name <- "ggplot-jitter-boxplot_total-expression_bulk-rna-types.jpg"
bulk.mean.rna.types.jpg.path <- here(save.path, bulk.mean.rna.types.jpg.name)
bulk.variance.rna.types.jpg.path <- here(save.path, bulk.variance.rna.types.jpg.name)
bulk.total.rna.types.jpg.path <- here(save.path, bulk.total.rna.types.jpg.name)
