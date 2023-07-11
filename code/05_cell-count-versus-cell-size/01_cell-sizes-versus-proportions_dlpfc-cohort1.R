
#
# Testing relationship between cell sizes and total cells.
#

# aggregation basis code:
# do.call(rbind, lapply(sample.id.vector, function(sample.id){}))

source("deconvo_method-paper/code/05_cell-count-versus-cell-size/00_parameters.R")
sapply(libv, library, character.only = T)

# load mae 
mae.filepath <- here("deconvo_method-paper", "outputs", "01_prepare-datasets", "mae_final.rda")
mae <- get(load(mae.filepath))
# get mae metadata for subsetting
mae.cd <- colData(mae.final)
sample.id.vector <- mae.cd$sample.id[complete.cases(mae.final)]

# get cell size tables
df.cellsize <- do.call(rbind, lapply(sample.id.vector, function(sample.id){
  mae.iter <- mae[,mae.cd$sample.id==sample.id,]
  table(mae.iter[3])
}))


