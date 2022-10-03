library(DeconvoBuddies)

# load example data
scef.fname <- "scef_kg-sg200_dlpfc_tran2021.rda"
scef <- get(load(scef.fname))

type.varname <- "cellType"

ngenes.byk <- 20

# test markerv
zsource_type(scef, type.varname)
# get markerdata
markerdata <- zsource_markerdata(zsource = scef, type.varname = type.varname)
dim(markerdata) # [1] 9188    8

