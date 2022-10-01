# 
# Contains parameters used by scripts in this task
#

# manage high-level paths for project
project.basepath <- "deconvo_method-paper"
code.dpath <- file.path(project.basepath, "code")

# project params
celltype.varname <- "cellType_broad_hc"
snrnaseq.fpath <- ""
bulkdata.fpath <- NA

# manage task-level paths
z.genes <- "filtered.markers"
param.paths <- paste0(task.id, "-00_param")
task.id <- "01"
task.dname <- "test-sk"

# task params
y.type <- "pseudobulk"
z.transform <- "sk"
z.postprocess <- NA
