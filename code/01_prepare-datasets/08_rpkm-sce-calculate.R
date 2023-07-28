library(edgeR)

# load data
sce <- mae[["sn1.rnaseq"]]
rse <- mae[["bulk.rnaseq"]]; rse <- rse[rownames(sce),]
sce.counts <- assays(sce)[["counts"]]
# load library sizes from complete snrnaseq object
sce.libsize.path <- "./deconvo_method-paper/outputs/01_prepare-datasets/total-counts-library-sizes_sce-dlpfc-cohort1.rda"
sce.libsize <- get(load(sce.libsize.path))
# set params
marker.lengths <- rowData(rse)$Length
library.sizes <- rep(1000, ncol(sce.counts))

# try with full sce object
sce.input <- DGEList(counts = sce.counts, 
                     lib.sizes = as.numeric(sce.libsize),
                     genes = data.frame(Length=rowData(rse)$Length))
#sce.input <- calcNormFactors(sce.input)
sce.rpkm <- rpkm(y = sce.input, log = F)
# returns error
#Error in if (median(f75) < 1e-20) { : 
#    missing value where TRUE/FALSE needed

#-----------------------------------
# try with z reference dataset -- k2
#-----------------------------------
# library size by k2 type
libsize.vector <- sapply(unique(sce[["k2"]]), function(cell.type){
  celltype.filter <- sce[["k2"]]==cell.type
  mean(sce.libsize[celltype.filter])})
sce.z.counts <- lute:::.get_z_from_sce(sce, "counts", "k2")
rpkm.input <- DGEList(counts = sce.z.counts, lib.sizes = libsize.vector,
                      genes = data.frame(Length=marker.lengths))
sce.z.rpkm <- rpkm(y = rpkm.input, log = F)
# save
sce.z.rpkm.path <- file.path("./deconvo_method-paper/outputs/01_prepare-datasets/sce-z-rpkm-k2_dlpfc-cohort1.rda")
save(sce.z.rpkm, file = sce.z.rpkm.path)

# check bulk rpkm
# assays(mae[["bulk.rpkm.rnaseq"]])[[1]]

#-----------------------------------
# try with z reference dataset -- k3
#-----------------------------------

#-----------------------------------
# try with z reference dataset -- k4
#-----------------------------------
