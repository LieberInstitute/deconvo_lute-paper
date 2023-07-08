#!/usr/bin/env R

libv <- c("SummarizedExperiment", "SingleCellExperiment")
sapply(libv, library, character.only = TRUE)

# manage paths
# save dir path
save.dname <- "deconvo_method-paper/outputs/"
save.dpath <- file.path(save.dname)
# sce rda path
sce.fname <- "sce_DLPFC.Rdata"
sce.dpath <- "DLPFC_snRNAseq/processed-data/sce"
sce.fpath <- file.path(sce.dpath, sce.fname)

# load data
sce <- get(load(sce.fpath))
pdat <- colData(sce)
pdat <- as.data.frame(pdat)
colnames(pdat)
# [1] "Sample"                "Barcode"               "key"
# [4] "SAMPLE_ID"             "pos"                   "BrNum"
# [7] "round"                 "Position"              "age"
# [10] "sex"                   "diagnosis"             "sum"
# [13] "detected"              "subsets_Mito_sum"      "subsets_Mito_detected"
# [16] "subsets_Mito_percent"  "total"                 "high_mito"
# [19] "low_sum"               "low_detected"          "discard_auto"
# [22] "doubletScore"          "prelimCluster"         "collapsedCluster"
# [25] "kmeans"                "sizeFactor"            "cellType_broad_k"
# [28] "cellType_k"            "cellType_broad_hc"     "cellType_hc"
# [31] "cellType_layer"        "layer_annotation"      "cellType_azimuth"

# save pdata
pdat.fname <- "pdata-sce_ro1-dlpfc.rda"
save(pdat, file = file.path(save.dpath, pdat.fname))

# get pdata by donor
pdat <- pdat[!duplicated(pdat$BrNum),]

# inspect sex by donor/brnum
table(pdat$sex, pdat$BrNum)
#     Br2720 Br2743 Br3942 Br6423 Br6432 Br6471 Br6522 Br8325 Br8492 Br8667
# F      1      0      0      0      0      0      0      1      1      1
# M      0      1      1      1      1      1      1      0      0      0


#-----------------
# randomize on sex
#-----------------
set.seed(0)

# get donors by sex
female.donorv <- unique(pdat[pdat$sex=="F",]$BrNum)
male.donorv <- unique(pdat[!pdat$BrNum %in% female.donorv,]$BrNum)

# get validation subset
val.f <- sample(female.donorv, 1)
val.m <- sample(male.donorv, 2)
valv <- c(val.f, val.m)
valv
# [1] "Br8667" "Br6522" "Br6432"

# get discovery subset
trainv <- unique(pdat[!pdat$BrNum %in% valv,]$BrNum)
trainv
# [1] "Br2720" "Br6471" "Br8492" "Br2743" "Br3942" "Br6423" "Br8325"
