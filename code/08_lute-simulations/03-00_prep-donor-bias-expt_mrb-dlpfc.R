#!/usr/bin/env R

# Author: Sean Maden
#
# Prepare datasets for donor bias experiments using the multi-region brain 
# DLPFC dataset.
#
#

libv <- c("scuttle", "glmGamPoi", "sva",
          "SingleCellExperiment", "SummarizedExperiment")
sapply(libv, library, character.only = TRUE)


#----------
# load data
#----------
# get save dpath
code.dname <- "08_lute-simulations"
proj.dname <- "deconvo_method-paper"
save.dpath <- file.path(proj.dname, "outputs", code.dname)
# get marker data
dfm.fname <- "markers-k2_db-mr2_sce-dlpfc-mrb.rda"
dfm <- get(load(file.path(save.dpath, dfm.fname)))
# get sce
sce.fname <- "sce-mrb_dlpfc.rda"
sce.fpath <- file.path(save.dpath, sce.fname)
sce <- get(load(sce.fpath))

#-------------------
# get top 20 markers
#-------------------
nmarker <- 20
typev <- unique(dfm$cellType.target)
dfmt <- do.call(rbind, lapply(typev, function(typei){
  dfm %>% filter(cellType.target == typei) %>% 
    arrange(rank_ratio) %>%
    top_n(n = nmarker)
}))

#------------------
# format assay data
#------------------
for(ii in names(assays(sce))){
  assays(sce)[[ii]] <- as.matrix(assays(sce)[[ii]])}

#---------------------
# perform downsampling
#---------------------
assays(sce)[["counts_ds"]] <- mexpr_downsample(assays(sce)[["counts"]])
assays(sce)[["logcounts_ds"]] <- mexpr_downsample(assays(sce)[["logcounts"]])

# save
sce.fname <- "sce_ds-assays_mrb-dlpfc.rda"
save(sce, file = file.path(save.dpath, sce.fname))

#---------------
# perform combat
#---------------
for(ai in names(assays(sce))){
  message("Working on assay: ", ai)
  ai <- "counts_ds"
  mi <- assays(sce)[[ai]]
  cnv <- colnames(mexpr)
  pheno <- data.frame(donor = sce[["donor"]], type = sce[["k2"]])
  mod <- model.matrix(~type, data = pheno)
  batch <- pheno$donor
  mi.adj <- ComBat(dat = mi, batch = batch, mod = mod, 
             par.prior = TRUE, prior.plots = FALSE)
  assays(sce)[[paste0(ai, "_combat")]] <- mi.adj
}

# save
sce.fname <- "sce_combat-ds-assays_mrb-dlpfc.rda"
save(sce, file = file.path(save.dpath, sce.fname))

#-------------------
# analysis functions
#------------------
get_anova_df <- function(mi, type = "complex"){
  lmi <- lm()
  avi <- anova(lmi)
  dfa <- as.data.frame(avi)
}

#-----------------------------
# get neg. binom. coefficients
#-----------------------------
# summarize by donor

# fit models

# get model coeffs

# save