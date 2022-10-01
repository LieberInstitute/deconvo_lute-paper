#!/usr/bin/env R

#
# Parameters for the project
#

# manage high-level paths for project
basepath.project <- "deconvo_method-paper"
code.dpath.project <- file.path(basepath.project, "code")
# project params
ktype.varname <- "cellType_broad_hc" # should be a valid variable name
zsource.fpath.project <- ""
ysource.fpath.project <- NA