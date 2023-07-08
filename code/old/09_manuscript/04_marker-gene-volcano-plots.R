# load full data
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
# get volcano plot
library(ggplot2); library(ggrepel)
filepath=""
dfp <- get(load(filepath))
dfp$y.axis <- -1*log10(dfp$p.value)
dfp$x.axis <- dfp$logFC
dfp.subset <- dfp[dfp$p.value <= 1e-3,]
newplot <- ggplot(dfp, aes(x = x.axis, y = y.axis)) + geom_point(alpha = 0.2) + 
geom_point(data = dfp.subset, aes(x = x.axis, y = y.axis), color = "red") +
geom_hline(yintercept = max(dfp.subset$y.axis), color = "black") +
geom_vline(xintercept = 0, color = "black") +
geom_label_repel(box.padding   = 0.35, point.padding = 0.8, segment.color = 'grey50')
# save plot object
ggsave(newplot)