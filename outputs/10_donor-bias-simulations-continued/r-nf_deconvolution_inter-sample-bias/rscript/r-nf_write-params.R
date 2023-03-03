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
                    default="/data/workflow-table.csv",
                    help = "Path to the workflow .csv table.")
parser$add_argument("-s", "--start_index", type="character", default="NULL",
                    help = "Index of start row to read from workflow table.")
parser$add_argument("-e", "--end_index", type="character", default="NULL",
                    help = "Index of final row to read from workflow table.")
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
start.index <- args$start_index
end.index <- args$end_index

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
  con <- file(pc.fpath)
  lp <- readLines(con)
} else{
  stop("Error, didn't find parameters .config table at: ", pc.fpath)
}

#---------------------------------
# parse workflow table row indices
#---------------------------------
max.row <- nrow(wt)
end.index <- ifelse(end.index > max.row, max.row, end.index)
cond <- is(start.index, "NULL")|is(end.index, "NULL")
cond <- cond|start.index=="NULL"|end.index=="NULL"
cond <- !cond & (start.index <= max.row)
cond <- cond & (start.index <= end.index)
if(cond){wt <- wt[start.index:end.index,]}
num.runs <- nrow(wt)
message("After parsing start and end row indices, found ", num.runs, " runs.")

#--------------------
# get new param lines
#--------------------
# get data to write
which.param <- which(grepl("param", md$type))
variable.names <- md[which.param,]$variable
wt.colnames <- colnames(wt)
variables.provided <- intersect(variable.names, wt.colnames)
message("In workflow table at ",wt.fpath,", found ", length(variables.provided), 
        " columns with annotations in params-metadata.csv.")

# make new lines
collapse.str <- ","
variable.lines <- lapply(variables.provided, function(variable.iter){
  message("working on provided variable '", variable.iter, "' ...")
  append.str <- md[variable.iter,]$append_string
  append.str <- ifelse(is.na(append.str)|append.str=="NA", "", append.str)
  values <- paste0(append.str, wt[,variable.iter])
  values <- paste0('"', gsub('"', '', values), '"', collapse = collapse.str)
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
message("Successfully wrote data for ", num.runs, 
        " runs to params file at ", pc.fpath, ".")
