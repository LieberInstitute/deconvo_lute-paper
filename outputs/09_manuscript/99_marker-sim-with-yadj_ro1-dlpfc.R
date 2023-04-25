library(lute)
library(dplyr)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(ggplot2)

sce <- lscef[[1]]
group.variable <- "Sample"
group.vector <- unique(sce[[group.variable]])

sample.main <- 2
sample.pb <- 4

sce.main <- sce[,sce[[group.variable]]==group.vector[sample.main]]
sce.pb <- sce[,sce[[group.variable]]==group.vector[sample.pb]]

z.main <- signature_matrix_from_sce(sce = sce.main, 
                                    celltype.variable = "k2", 
                                    summary.method = "mean", 
                                    assay.name = "counts")
zpb <- signature_matrix_from_sce(sce = sce.pb, 
                                    celltype.variable = "k2", 
                                    summary.method = "mean", 
                                    assay.name = "counts")
p <- c(0.2, 0.8)
s <- c(3, 10)
zs <- sweep(zpb, 2, s, FUN = "*")
ypb <- t(t(p) %*% t(zs))
ypb <- ypb + 10
rownames(ypb) <- rownames(sce)

# without adj
# result1 <- lute(z = z.main, y = ypb)
result1 <- deconvolution(nnlsParam(z = z.main, y = ypb, s = c(1,1)))
result1/sum(result1)

# with adj
#ypb.main <- t(t(p) %*% t(z.main))
#data.table <- cbind(ypb, ypb.main) %>% as.data.frame()
#colnames(data.table) <- c("ypb", "ymain")
data.table <- cbind(ypb, z.main) %>% as.data.frame()
colnames(data.table) <- c("ypb", "z1", "z2")
lm <- lm(ypb ~ z1 + z2, data = data.table)
y.fit <- matrix(lm$fitted.values, ncol = 1)
rownames(y.fit) <- rownames(sce)
result2 <- deconvolution(nnlsParam(z = z.main, y = y.fit, s = c(1,1)))
result2/sum(result2)

result3 <- deconvolution(nnlsParam(z = z.main, y = y.fit, s = c(3, 10)))
result3/sum(result3)

result4 <- deconvolution(nnlsParam(z = z.main, y = ypb, s = c(3, 10)))
result4/sum(result4)

nnls::nnls(z.main, ypb)

ggplot(plot.data, aes(x = glial, y = neuron)) + geom_point() + geom_abline(slope = 1)


