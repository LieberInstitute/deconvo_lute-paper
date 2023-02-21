#!/usr/bin/env R

# url <- "https://github.com/cozygene/bisque"
# devtools::install_github(url)
# library(BisqueRNA)

install.packages("BisqueRNA")
library(BisqueRNA)

# parse sce data
subject.variable <- "Sample"
# parse sce metadata
cd <- colData(sce)
cdf <- as.data.frame(cd)
# adf <- annotatedDataFrameFrom(list(cdf), byrow=FALSE)
sc.mexpr <- as.matrix(assays(sce)[[assay.name]])
unique.subjects <- unique(sce[[subject.variable]])
if(length(unique.subjects)==1){
  # append expr
  sc.mexpr1 <- sc.mexpr
  colnames(sc.mexpr1) <- paste0(colnames(sc.mexpr1), "_rep")
  sc.mexpr <- cbind(sc.mexpr, sc.mexpr1)
  # append coldata
  cdf[,subject.variable] <- "1"
  cdf1 <- cdf
  rownames(cdf1) <- colnames(sc.mexpr1)
  cdf1[,subject.variable] <- "2"
  cdf <- rbind(cdf, cdf1)
  # parse bulk data
  bulk.mexpr <- cbind(ypb, ypb, ypb)
  colnames(bulk.mexpr) <- c(unique.subjects, "1", "2")
}
adf <- AnnotatedDataFrame(cdf)
sc.eset <- ExpressionSet(assayData = sc.mexpr, phenoData = adf)
bulk.eset <- ExpressionSet(assayData = bulk.mexpr)
res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset,
                                              cell.types = celltype.variable,
                                              subject.names = subject.variable)
