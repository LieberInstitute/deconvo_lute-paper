#!/usr/bin/env R

#
# Main script to parse parameters for a task. For a given task, this reads in and
# parses the project and task parameters at expected locations (e.g. ./code/* for 
# project parameters file, and ./code/taskdir/* for task parameters file).
# 
# Notes:
# * Recall we use standard terminology: Y = pi * Z
#

param.proj <- file.path("deconvo_method-paper", "code", "00_parameters.R")
param.task <- file.path("deconvo_method-paper", "code", "01_test-sk")
source.dpath <- file.path("deconvo_method-paper", "source")

# get params
source(proj.params.fpath)
source(task.params.fpath)

# source functions, etc. shared across tasks
for(fn in list.files(source.dpath)){source(source.dpath)}


# parse pi options
pi.data <- "NA"
if(is.na(pi.source.task)){
	if(y.type.task=="pseudobulk"){
		message("Using pseudobulk for pi.")
		pi.data <- se.pb$pidata
	}
} else{
	pi.data <- get(load(pi.source.fpath.task))
}

# parse z options
# parse z data options
z.data <- "NA"
message("parsing options for z source data...")
if(is.na(z.data.fpath.task)){
	message("skipping z source for task...")
} else{
	z.data <- get(load(z.data.fpath.task))
}
# parse z source options
z.source <- "NA"
if(!is.na(z.source.fpath.task)){
	message("loading sce object...")
	z.source <- get(load(z.source.fpath.task))
} else{
	message("skipping z source...")
}








# parse y options
message("parsing options for y mixed signals matrix...")
if(!is.na(ysource.fpath.project)){
	y <- get(load(ysource.fpath.project))
	} else if(y.type.task == "pseudobulk"){
		message("making pseudobulk dataset")
	} else{
		message("no y data to parse.")
	}
