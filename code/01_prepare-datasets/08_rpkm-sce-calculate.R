library(edgeR)


sce <- mae[["sn1.rnaseq"]]
rse <- mae[["bulk.rnaseq"]]; rse <- rse[rownames(sce),]

y <- DGEList(counts=assays(sce)[["counts"]],genes=data.frame(Length=rowData(rse)$Length))
y <- calcNormFactors(y)

RPKM <- rpkm(y)
sce.rpkm <- rpkm(assays(sce1)[["counts"]])
