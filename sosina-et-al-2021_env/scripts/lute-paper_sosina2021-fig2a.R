#!/usr/bin/env R

# Author: Sean Maden
#
# Code to replicate Figure 2a from Sosina et al 2021 (DOI:10.12688/f1000research.50858.1)
#
#


scripts.path <- "scripts"
libv <- c("here", "SingleCellExperiment", "data.table", "dplyr", "lute", "ggplot2", "MuSiC")
sapply(libv, library, character.only = T)

#----------
# load data
#----------
# load bulk expression
rse.path <- file.path("files", "NAc_rse_gene_withCompEsts_update.rda")
rse.bulk <- get(load(rse.path))
# use gene symbols
rownames(rse.bulk) <- rowData(rse.bulk)$Symbol

# load proportions
prop.path <- file.path("files", "proportions_sm.csv")
prop <- fread(prop.path, data.table = F)
prop <- prop[,c(1:2)]

#----------------------
# read in counts matrix
#----------------------
path <- file.path("files", "countMatrix_n4169-NAc-nuclei.csv")
counts <- fread(path, data.table = F)
rownames(counts) <- as.character(counts[,1])
counts <- counts[,c(2:ncol(counts))]

#--------------------------
# get reference from counts
#--------------------------
pheno.path <- here("files", "pd-cellTypeAssignment_n4169.csv")
pheno <- fread(pheno.path, data.table = F)
identical(colnames(counts), rownames(pheno))

# check
nrow(pheno) == ncol(counts) # TRUE
head(pheno[,1])
head(colnames(counts))

pheno.reorder <- order(match(pheno[,1], colnames(counts)))
pheno <- pheno[pheno.reorder,]
identical(as.character(pheno[,1]), as.character(colnames(counts)))

# save sce of nac 2 donors
sce.nac <- SingleCellExperiment(assays = list(counts = counts), colData = pheno)
save(sce.nac, file = "sce_nac-2donors-sosina-strategies_lieber.rda")

# get zref using means, then save
zref <- lute:::.get_z_from_sce(sce = sce.nac, "counts", "nucleusCellType")
save(zref, file = "zref_nac-2donors-sosina-strategies_lieber.rda")

# load markers
top50.path <- here("files", "clusterMarkers_50-per-cluster_n4169.csv")
top50 <- fread(top50.path, data.table = F)
top25.path <- here("files", "clusterMarkers_25-per-cluster_n4169.csv")
top25 <- fread(top25.path, data.table = F)

# get marker subsets for y, z
# markers in both y,z
top50.yz <- intersect(top50$gene, intersect(rownames(rse.bulk), rownames(zref)))
top25.yz <- intersect(top25$gene, intersect(rownames(rse.bulk), rownames(zref)))
y50 <- rse.bulk[top50.yz,]
y25 <- rse.bulk[top25.yz,]
z50 <- zref[top50.yz,]
z25 <- zref[top25.yz,]

#------------------
# music predictions
#------------------
# prep
sc.sce.music <- sce.nac[rownames(z50),]
sc.sce.music$Neuron <- ifelse(sc.sce.music$nucleusCellType == "Neuron", 
                              TRUE, FALSE)
sc.sce.music$Sample <- gsub("\\..*", "", sc.sce.music$V1)

# predict
music_est <- music_prop(bulk.mtx = assays(y50)[["counts"]], 
                        sc.sce = sc.sce.music, 
                        clusters = "Neuron",
                        samples = "Sample")

# plot
pred <- music_est$Est.prop.weighted  %>% as.data.frame()
identical(rownames(pred), prop$samples)
prop$music <- pred[,"TRUE"]
dfp <- data.frame(true = prop[,2], music = pred$`TRUE`)

ggpt.fig2a <- ggplot(dfp, aes(x = music, y = true)) + 
  geom_point(alpha = 1, shape = "square", color = "forestgreen", size = 3) +
  geom_abline(slope = 1, intercept = 0) +
  xlim(0, 0.5) + ylim(0, 0.4) + theme_classic()

jpeg("lute-replicate-sosina2021-fig2a.jpg", width = 5, height = 5, units = "in", res = 400)
ggpt.fig2a
dev.off()

#-----------------
# save environment
#-----------------
save.image(file = file.path(env.path, "env_predict-music.RData"))
