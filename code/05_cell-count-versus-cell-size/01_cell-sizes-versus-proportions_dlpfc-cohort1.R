
#
# Testing relationship between cell sizes and total cells.
#

# load mae 
mae.filepath <- here("deconvo_method-paper", "outputs", "01_prepare-datasets", "mae_final.rda")
load(mae.filepath)

# get cell size tables
mae.final[1,1,1]
