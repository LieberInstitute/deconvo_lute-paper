library(edgeR)

# load data
sce <- mae[["sn1.rnaseq"]]
rse <- mae[["bulk.rnaseq"]]; rse <- rse[rownames(sce),]
sce.counts <- assays(sce)[["counts"]]

# set params
marker.lengths <- rowData(rse)$Length
library.sizes <- rep(1000, ncol(sce.counts))

# try with full sce object
sce.input <- DGEList(counts = sce.counts, 
                     lib.sizes = rep(1000, ncol(sce.counts)),
                     genes = data.frame(Length=rowData(rse)$Length))
sce.input <- calcNormFactors(sce.input)
sce.rpkm <- rpkm(y = sce.input)
# returns error
#Error in if (median(f75) < 1e-20) { : 
#    missing value where TRUE/FALSE needed

# try with reference dataset
sce.z.counts <- lute:::.get_z_from_sce(sce, "counts", "k2")
rpkm.input <- DGEList(counts = sce.z.counts, 
                     lib.sizes = library.sizes,
                     genes = data.frame(Length=marker.lengths))
rpkm.input <- calcNormFactors(sce.input)
sce.z.rpkm <- rpkm(y = rpkm.input)
# returns error
# Error in if (median(f75) < 1e-20) { : 
#    missing value where TRUE/FALSE needed
