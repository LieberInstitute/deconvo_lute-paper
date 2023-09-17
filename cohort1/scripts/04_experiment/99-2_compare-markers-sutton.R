#!/usr/bin/env R

# Compare marker gene set overlaps.
#
#
#


# sigsBrain.rda source: https://github.com/Voineagulab/DeconRNAShiny/blob/main/sigsBrain.rda
sutton.markers.path <- file.path("deconvo_method-paper/outputs/15_k2-simulations_within-sample-matched/sigsBrain.rda")
sutton.markers <- get(load(sutton.markers.path))

gene.info.path <- file.path("deconvo_method-paper/outputs/15_k2-simulations_within-sample-matched/geneInfo.rda")
gene.info <- get(load(gene.info.path))

# map ensg to symbol
org.Hs.eg.db.

library("AnnotationDbi")
library("org.Hs.eg.db")
res <- results(dds,alpha=.05, contrast=c("Type", "Disease", "Control"))
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="SYMBOL", 
                     keytype="ENSEMBL", 
                     multiVals="first")

library(biomaRt)
mart<- useDataset("hsapiens_gene_ensembl", useMart("ENSEMBL_MART_ENSEMBL"))

dim(sigsBrain$MB) # [1] 8580    5
ensemblGenes <- rownames(sigsBrain$MB)
ens <- getBM(filters="ensembl_gene_id", 
             attributes=c("ensembl_gene_id", "hgnc_symbol"),
             values=ensemblGenes,
             mart=mart)
length(intersect(ens$hgnc_symbol, markers.vector.overall.top)) # 34, 31

dim(sigsBrain$LK) # [1] 9577    8
ensemblGenes <- rownames(sigsBrain$LK)
ens <- getBM(filters="ensembl_gene_id", 
             attributes=c("ensembl_gene_id", "hgnc_symbol"),
             values=ensemblGenes,
             mart=mart)
length(intersect(ens$hgnc_symbol, markers.vector.overall.top)) # 34
# 
data <- sigsBrain$LK
glia.vector <- c("Astrocytes", "Microglia", "Oligodendrocytes", "OPCs")
data$glia <- rowMedians(as.matrix(data[,glia.vector])) + 1
data$mean.ratios <- data$Neurons/data$glia
data <- data[rev(order(data$mean.ratios)),]
ensemblGenes <- rownames(data)[seq(1000)]
ens <- getBM(filters="ensembl_gene_id", 
             attributes=c("ensembl_gene_id", "hgnc_symbol"),
             values=ensemblGenes,
             mart=mart)
length(intersect(ens$hgnc_symbol, markers.vector.overall.top)) # 1
intersect(ens$hgnc_symbol, markers.vector.overall.top) # [1] "PODXL2"

# TS
data <- sigsBrain$TS
dim(sigsBrain$TS) # [1] 11582     8
ensemblGenes <- rownames(sigsBrain$TS)
ens <- getBM(filters="ensembl_gene_id", 
             attributes=c("ensembl_gene_id", "hgnc_symbol"),
             values=ensemblGenes,
             mart=mart)
length(intersect(ens$hgnc_symbol, markers.vector.overall.top)) # 33
# top 100 neurons
data <- sigsBrain$TS
glia.vector <- c("Astrocytes", "Microglia", "Oligodendrocytes", "OPCs")
data$glia <- rowMedians(as.matrix(data[,glia.vector])) + 1
data$mean.ratios <- data$Neurons/data$glia
data <- data[rev(order(data$mean.ratios)),]
ensemblGenes <- rownames(data)[seq(1000)]
ens <- getBM(filters="ensembl_gene_id", 
             attributes=c("ensembl_gene_id", "hgnc_symbol"),
             values=ensemblGenes,
             mart=mart)
length(intersect(ens$hgnc_symbol, markers.vector.overall.top)) # 3
intersect(ens$hgnc_symbol, markers.vector.overall.top) # "RASGRF1" "PODXL2"  "CA10"

# ca
# top 100 mean ratios
data <- sigsBrain$CA
glia.vector <- c("Astrocytes", "Microglia", "Oligodendrocytes", "OPCs")
data$glia <- rowMedians(as.matrix(data[,glia.vector])) + 1
data$mean.ratios <- data$Neurons/data$glia
data <- data[rev(order(data$mean.ratios)),]
ensemblGenes <- rownames(data)[seq(1000)]
ens <- getBM(filters="ensembl_gene_id", 
             attributes=c("ensembl_gene_id", "hgnc_symbol"),
             values=ensemblGenes,
             mart=mart)
length(intersect(ens$hgnc_symbol, markers.vector.overall.top)) # 3
intersect(ens$hgnc_symbol, markers.vector.overall.top) # "RASGRF1" "CA10"    "RANBP17"

# NG
data <- sigsBrain$NG
glia.vector <- c("Astrocytes", "Microglia", "Oligodendrocytes", "OPCs")
data$glia <- rowMedians(as.matrix(data[,glia.vector])) + 1
data$mean.ratios <- data$Neurons/data$glia
data <- data[rev(order(data$mean.ratios)),]
ensemblGenes <- rownames(data)[seq(1000)]
ens <- getBM(filters="ensembl_gene_id", 
             attributes=c("ensembl_gene_id", "hgnc_symbol"),
             values=ensemblGenes,
             mart=mart)
length(intersect(ens$hgnc_symbol, markers.vector.overall.top)) # 3
intersect(ens$hgnc_symbol, markers.vector.overall.top) # "RASGRF1" "MTFMT"

