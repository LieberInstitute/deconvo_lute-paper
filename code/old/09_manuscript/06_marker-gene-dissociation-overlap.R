# load full expression data
sce.fname <- "sce_DLPFC.Rdata"
sce.dpath <- "DLPFC_snRNAseq/processed-data/sce"
sce.fpath <- file.path(sce.dpath, sce.fname)
sce <- get(load(sce.fpath))
# load marker data
markers.filename <- "list-scef_markers-k2-k3-k4_ro1-dlpfc.rda"
sce.markers.path <- file.path("deconvo_method-paper", "outputs", 
	"09_manuscript", markers.filename)
lsce <- get(load(sce.markers.path))
sce <- lsce[["k2"]]
rm(lsce)
# load dissociation markers
# get overlaps with markers
# get overlaps with full
# compare fractions overlapping
# 