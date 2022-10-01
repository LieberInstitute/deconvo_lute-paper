#!/usr/bin/env R

#
# Main script to parse parameters for a task
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




# parse z options
message("parsing options for z source data...")
if(is.na(zsource.fpath.task)){
	message("skipping z source for task...")
	
}



if(!is.na(zsource.fpath.project)){
	message("loading sce object...")
	sce <- get(load(zsource.fpath.project))
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
