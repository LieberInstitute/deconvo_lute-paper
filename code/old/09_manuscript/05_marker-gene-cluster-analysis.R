# load marker data
markers.filename <- "list-scef_markers-k2-k3-k4_ro1-dlpfc.rda"
sce.markers.path <- file.path("deconvo_method-paper", "outputs", 
	"09_manuscript", markers.filename)
lsce <- get(load(sce.markers.path))
sce <- lsce[["k2"]]
rm(lsce)
# cluster cells
# get plot data
# plot cluster data
