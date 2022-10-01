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
# pu-specific parameters
pi.source.task <- NA
#
#
#
# y-specific parameters
y.type.task <- "pseudobulk"
#
#
#
#
#
# z-specific parameters
z.fpath.task <- ""
zsource.fpath.task <- NA
z.transform.task <- "sk"
z.postprocess.task <- NA
z.type.task <- "mean.mr.genemarkers"