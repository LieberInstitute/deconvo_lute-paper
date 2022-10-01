#!/usr/bin/env R

# 
# Contains parameters used by scripts in this task
#

# manage task-level paths
z.genes.task <- "filtered.markers"
param.paths.task <- paste0(task.id, "-00_param")
id.task <- "01"
dname.task <- "test-sk"

# task params
#
#
#
# pi parameters
# note: if this is NA and y.type.task == "pseudobulk", takes pi from pseudobulk data.
pi.source.task <- NA 
pi.source.fpath.task <- ""
#
#
#
# y parameters
y.type.task <- "pseudobulk"
#
#
#
#
#
# z parameters
# source data
# note: if na, ignores source data
# note: options here, either NA or z.source.fpath.project
z.source.fpath.task <- NA
# z
# note: if na, ignores z
z.data.fpath.task <- "path_to_z_object" 
# z.type.task <- "mean.mr.genemarkers"
# transformations
# note: parameters here will source transformation functions at `./source...`
z.transform.task <- "sk"
z.postprocess.task <- NA
