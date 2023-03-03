#!/usr/bin/env R

# url <- "https://github.com/cozygene/bisque"
# devtools::install_github(url)
# library(BisqueRNA)

install.packages("BisqueRNA")
library(BisqueRNA)


celltype.variable <- "k2"
sample.variable <- "Sample"
assay.name <- "counts_adj"

# for between-sample data
unique.samples <- unique(sce[[sample.variable]])[1:3]
sce <- sce[,sce[[sample.variable]] %in% unique.samples]
# bind new ypb column for each sample
unique.samples <- unique(sce[[sample.variable]])
num.col.bulk <- length(unique.samples) + 1
bulk.mexpr <- do.call(cbind, lapply(seq(num.col.bulk), function(ii){ypb}))
colnames(bulk.mexpr) <- c(unique.samples, "rep")
bulk.eset <- ExpressionSet(assayData = bulk.mexpr)
# parse sce data
cd <- colData(sce)
cdf <- cdf2 <- as.data.frame(cd)
adf <- AnnotatedDataFrame(cdf)
sc.mexpr <- as.matrix(assays(sce)[[assay.name]])
if(length(unique.samples)==1){
  # force order
  var.lvl <- unique(cdf[,celltype.variable])
  cdf[,celltype.variable] <- factor(cdf[,celltype.variable], 
                                    levels=var.lvl[order(var.lvl)])
  cdf2[,sample.variable] <- "rep"
  rownames(cdf2) <- paste0(rownames(cdf2), "_2")
  cdf <- rbind(cdf, cdf2)
  adf <- AnnotatedDataFrame(cdf)
  sc.mexpr <- sc.mexpr2 <- as.matrix(assays(sce)[[assay.name]])
  colnames(sc.mexpr2) <- paste0(colnames(sc.mexpr2), "_2")
  sc.mexpr <- cbind(sc.mexpr, sc.mexpr2)
}
# filter 0 expression
maxv <- colMaxs(sc.mexpr)
sc.mexpr <- sc.mexpr[,!maxv==0]
adf <- adf[!maxv==0,]
sc.eset <- ExpressionSet(assayData = sc.mexpr, phenoData = adf)
# get proportions
res <- ReferenceBasedDecomposition(bulk.eset, sc.eset, 
                                   cell.types = celltype.variable,
                                   subject.names = sample.variable,
                                   use.overlap = FALSE)$bulk.props[,"rep"]

# for within-sample data
# append data to meet minimum sample requirements
# use overlaps
# use original data
cd <- colData(sce)
cdf <- cdf2 <- as.data.frame(cd)
# force order
var.lvl <- unique(cdf[,celltype.variable])
cdf[,celltype.variable] <- factor(cdf[,celltype.variable], 
                                  levels=var.lvl[order(var.lvl)])
cdf2[,sample.variable] <- "rep"
rownames(cdf2) <- paste0(rownames(cdf2), "_2")
cdf <- rbind(cdf, cdf2)
adf <- AnnotatedDataFrame(cdf)
sc.mexpr <- sc.mexpr2 <- as.matrix(assays(sce)[[assay.name]])
colnames(sc.mexpr2) <- paste0(colnames(sc.mexpr2), "_2")
sc.mexpr <- cbind(sc.mexpr, sc.mexpr2)
sc.eset <- ExpressionSet(assayData = sc.mexpr, phenoData = adf)
bulk.mexpr <- cbind(ypb, ypb, ypb)
colnames(bulk.mexpr) <- c(unique(cdf[,sample.variable]), "2")
bulk.eset <- ExpressionSet(assayData = bulk.mexpr)
res <- ReferenceBasedDecomposition(bulk.eset, sc.eset, 
                                   cell.types = celltype.variable,
                                   subject.names = sample.variable,
                                   use.overlap = TRUE)$bulk.props[,"2"]

# append data to meet minimum sample requirements
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
res <- ReferenceBasedDecomposition(bulk.eset, sc.eset, 
                                   cell.types = celltype.variable,
                                   subject.names = subject.variable,
                                   use.overlap = FALSE)



