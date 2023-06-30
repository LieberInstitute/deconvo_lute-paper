#!/usr/bin/env R

# Author: Sean Maden
#
# Code to replicate Figure 2b from Sosina et al 2021 (DOI:10.12688/f1000research.50858.1)
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
# cell size factors
#------------------
# source(file.path(scripts.path, "chunk2_prepare-cellsize-osm.R"))
# scripts.path <- "scripts"
# source(file.path(scripts.path, "chunk_prepare-cellsizes-osm.R"))

cell_size_osm_fish_cellarea <- c("90.87", "122.96")
cell_size_osm_fish_nrna_intensity <- c("180.46", "198.86")

# make dataframes
osm.nrna <- data.frame(celltypes = c(FALSE, TRUE), 
                       cellsizes = as.numeric(cell_size_osm_fish_nrna_intensity))
osm.cellarea <- data.frame(celltypes = c(FALSE, TRUE), 
                           cellsizes = as.numeric(cell_size_osm_fish_cellarea))

# source(file.path(scripts.path,"chunk_predict-plot-music-cellsize-osm.R"))
#------------------
# music predictions
#------------------
# prep
sc.sce.music <- sce.nac[rownames(z50),]
sc.sce.music$Neuron <- ifelse(sc.sce.music$nucleusCellType == "Neuron", 
                              TRUE, FALSE)
sc.sce.music$Sample <- gsub("\\..*", "", sc.sce.music$V1)

# predict
music_est_nrna <- music_prop(bulk.mtx = assays(y50)[["counts"]], 
                             sc.sce = sc.sce.music, 
                             clusters = "Neuron",
                             samples = "Sample",
                             cell_size = osm.nrna)

music_est_cellarea <- music_prop(bulk.mtx = assays(y50)[["counts"]], 
                                 sc.sce = sc.sce.music, 
                                 clusters = "Neuron",
                                 samples = "Sample",
                                 cell_size = osm.cellarea)

# plot
pred.nrna <- music_est_nrna$Est.prop.weighted  %>% as.data.frame()
pred.cellarea <- music_est_cellarea$Est.prop.weighted  %>% as.data.frame()

identical(rownames(pred.nrna), prop$samples)
identical(rownames(pred.cellarea), prop$samples)
prop$music.osm.nrna <- pred.nrna[,"TRUE"]
prop$music.osm.cellarea <- pred.cellarea[,"TRUE"]

dfp <- data.frame(true = prop[,2], 
                  music.osm.nrna = pred.nrna$`TRUE`, 
                  music.osm.cellarea = pred.cellarea$`TRUE`)

ggpt.fig2b.osm.nrna <- 
  ggplot(dfp, aes(x = music.osm.nrna, y = true)) + 
  geom_point(alpha = 1, col = "springgreen4", shape = 15) +
  geom_abline(slope = 1, intercept = 0) +
  xlim(0, 0.5) + ylim(0, 0.4) + theme_classic()

ggpt.fig2b.osm.cellarea <- 
  ggplot(dfp, aes(x = music.osm.cellarea, y = true)) + 
  geom_point(alpha = 1, col = "springgreen4", shape = 15) +
  geom_abline(slope = 1, intercept = 0) +
  xlim(0, 0.5) + ylim(0, 0.4) + theme_classic()

save.image(file = file.path(env.path, "env_predict-music.RData"))
