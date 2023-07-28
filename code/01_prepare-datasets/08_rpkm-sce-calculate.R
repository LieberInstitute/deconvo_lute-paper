library(edgeR)

# load data
sce <- mae[["sn1.rnaseq"]]
rse <- mae[["bulk.rnaseq"]]; rse <- rse[rownames(sce),]
sce.counts <- assays(sce)[["counts"]]
# load library sizes from complete snrnaseq object

# set params
marker.lengths <- rowData(rse)$Length
library.sizes <- rep(1000, ncol(sce.counts))

# try with full sce object
sce.input <- DGEList(counts = sce.counts, 
                     lib.sizes = colSums(sce.counts)+1,
                     genes = data.frame(Length=rowData(rse)$Length))
sce.input <- calcNormFactors(sce.input)
sce.rpkm <- rpkm(y = sce.input, log = F)
# returns error
#Error in if (median(f75) < 1e-20) { : 
#    missing value where TRUE/FALSE needed

# try with reference dataset
sce.z.counts <- lute:::.get_z_from_sce(sce, "counts", "k2")
rpkm.input <- DGEList(counts = sce.z.counts, 
                     lib.sizes = colSums(sce.z.counts),
                     genes = data.frame(Length=marker.lengths))
sce.z.rpkm <- rpkm(y = rpkm.input, log = F)
