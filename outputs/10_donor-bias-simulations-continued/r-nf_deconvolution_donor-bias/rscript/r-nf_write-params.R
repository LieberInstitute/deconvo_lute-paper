#!/usr/bin/env R

# Author: Sean Maden
#
# Write params.config from a CSV file
#

libv <- c("argparse")
sapply(libv, library, character.only = T)

#----------------
# parse arguments
#----------------
parser <- ArgumentParser() # create parser object
# data arguments
parser$add_argument("-f", "--workflow_table_filepath", type="character", 
                    default="workflow-table.csv",
                    help = "Path to the workflow .csv table.")
parser$add_argument("-m", "--parameters_metadata", type="character", 
                    default="params-metadata.csv",
                    help = "Path to the parameters metadata .csv table.")
parser$add_argument("-p", "--parameters_config", type="character", 
                    default="params.config",
                    help = "Path to the parameters .config file.")

args <- parser$parse_args() # get parser object

# parse provided arguments
wt.fpath <- args$workflow_table_filepath
pm.fpath <- args$parameters_metadata
pc.fpath <- args$parameters_config

#------------
# load data
#------------
# load table data
# workflow table
if(file.exists(wt.fpath)){
  wt <- read.csv(wt.fpath)
} else{
  stop("Error, didn't find workflow table at: ", wt.fpath)
}
# parameters metadata
if(file.exists(pm.fpath)){
  md <- read.csv(pm.fpath)
} else{
  stop("Error, didn't find parameters metadata table at: ", pm.fpath)
}
# parameters config
if(file.exists(pc.fpath)){
  lp <- readLines(pc.fpath)
} else{
  stop("Error, didn't find parameters .config table at: ", pc.fpath)
}

#--------------------
# get new param lines
#--------------------
# get data to write
which.param <- which(md$type=="param")
variable.names <- md[which.param,]$variable
wt.colnames <- colnames(wt)
variables.provided <- intersect(variable.names, wt.colnames)

# make new lines
variable.lines <- lapply(variables.provided, function(variable.iter){
  append.str <- md[variable.iter,]$append_string
  append.str <- ifelse(is.na(append.str)|append.str=="NA", "", append.str)
  values <- paste0(append.str, wt[,variable.iter])
  values <- paste0('"', gsub('"', '', values), '"', collapse = ", ")
  paste0("    ", variable.iter, " = [", values, "]")
})

#---------------------
# update params.config
#---------------------
# get start line
line.start.str <- ".*// LIST CHANNEL INPUTS"
which.start <- which(grepl(line.start.str, lp))

# write params
new.lines <- lp[1:(which.start+1)]
new.lines <- c(new.lines, unlist(variable.lines), "}")
writeLines(text = new.lines, con = con)
