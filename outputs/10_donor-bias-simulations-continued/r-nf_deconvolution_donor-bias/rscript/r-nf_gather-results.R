#!/usr/bin/env R

# Author: Sean Maden
#
# Gather deconvolution results from run.
#
#


libv <- c("argparse")
sapply(libv, library, character.only = T)

#----------------
# parse arguments
#----------------
parser <- ArgumentParser() # create parser object
# data arguments
parser$add_argument("-t", "--min_timestamp_filter", type="character", 
                    default="NULL",
                    help = paste0("Filter for latest timestamp. ",
                                  "If provided, only use results ",
                                  "after specified time. If NULL, use all ",
                                  "identified results"))
parser$add_argument("-s", "--string_identifier", type="character", 
                    default='deconvolution_analysis_.*',
                    help = paste0("Regex char string for identifying ",
                                  "results to gather."))
parser$add_argument("-d", "--results_directory", 
                    type="character", default='results',
                    help = paste0("Name of the publish directory ",
                                  "containing results"))
parser$add_argument("-m", "--deconvolution_method_colname", 
                    type="character", default='deconvolution_method',
                    help = paste0("Column name containing deconvolution method"))
parser$add_argument("-l", "--celltype_labels_colname", 
                    type="character", default='type_labels',
                    help = paste0("Column name containing celltype labels."))
parser$add_argument("-p", "--string_predicted_proportions", 
                    type="character", default='^prop\\.pred\\.',
                    help = paste0("Regex char string for predicted ",
                                  "cell type proportions column names."))
parser$add_argument("-r", "--string_true_proportions", 
                    type="character", default='^prop\\.true\\.',
                    help = paste0("Regex char string for true ",
                                  "cell type proportions column names."))
                    

args <- parser$parse_args() # get parser object

# parse provided arguments
timestamp <- args$min_timestamp_filter
identifier <- args$string_identifier
publish.dir <- args$results_directory
method.cname <- args$deconvolution_method_cname
type.labels.cname <- args$celltype_labels_cname
filt.str.pred <- args$string_predicted_proportions
filt.str.true <- args$string_true_proportions

#----------------
# parse filenames
#----------------
lfv <- list.files(publish.dir)
lfv.filt <- lfv[grepl(fname.str, lfv)]

#-------------------
# load and bind data
#-------------------
dfres <- do.call(rbind, lapply(lfv.filt, function(fni){
  message(fni)
  dfi <- read.csv(file.path(publish.dir, fni))
  dfi[,2:ncol(dfi)]
}))

#--------------------
# append rmse by type
#--------------------
res.colnames <- colnames(dfres)
if(type.labels.cname %in% colnames(dfres)){
  message("Getting RMSE within types...")
  cnvf <- colnames(dfres)
  cnvf <- cnvf[grepl(filt.str.pred, cnvf)]
  unique.types <- unique(gsub(".*\\.", "", cnvf))
  df.rmse <- do.call(cbind, lapply(unique.types, function(typei){
    cname.pred.str <- paste0(filt.str.pred, typei)
    cname.true.str <- paste0(filt.str.true, typei)
    pred.prop <- dfres[,grepl(cname.pred.str, res.colnames)][1]
    true.prop <- dfres[,grepl(cname.true.str, res.colnames)][1]
    pred.prop <- as.numeric(pred.prop)
    true.prop <- as.numeric(true.prop)
    rmsei <- sqrt(mean((pred.prop-true.prop)^2))
    rep(rmsei, nrow(dfres))
  }))
  colnames(df.rmse) <- paste0("rmse.", unique.types)
  dfres <- cbind(dfres, df.rmse)
} else{
  message("Didn't find any columns called '",type.labels.cname,"'. ",
          "Skipping within-type RMSE calculation...")
}

#----------------------
# append rmse by method
#----------------------
if(method.cname %in% colnames(dfres)){
  message("Getting RMSE by method...")
  cnvf <- colnames(dfres)
  cnvf <- cnvf[grepl(method.cname, cnvf)]
  unique.methods <- unique(dfres[,cnvf])
  df.rmse.method <- do.call(cbind, lapply(unique.methods, function(methodi){
    method.filt <- dfres[,method.cname]==methodi
    dff <- dfres[method.filt,]
    truev <- unlist(dff[,grepl(filt.str.true, colnames(dff))])
    predv <- unlist(dff[,grepl(filt.str.pred, colnames(dff))])
    truev <- as.numeric(truev)
    predv <- as.numeric(predv)
    rmsei <- sqrt(mean((predv-truev)^2))
    rep(rmsei, nrow(dfres))
  }))
  colnames(df.rmse.method) <- paste0("rmse.", unique.methods)
  dfres <- cbind(dfres, df.rmse.method)
} else{
  message("Didn't find any columns called '",method.cname,"'. ",
          "Skipping within-method RMSE calculation...")
}

#-------------
# save results
#-------------
# get results filename
ts <- as.character(as.numeric(Sys.time()))
new.filename <- paste0("results_table_",ts,".csv")

# save results table
write.csv(dfres, file = new.filename, row.names = FALSE)
