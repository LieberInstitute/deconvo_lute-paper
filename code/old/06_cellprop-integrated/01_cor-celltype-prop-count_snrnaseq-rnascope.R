#!/usr/bin/env R

# Author: Sean Maden
#
# Compare initial cell proportions from snRNA-seq and RNAscope.
#
# This script compares cell proportion, count correlations across the binned 
# data (i.e. median values), with the following details:
# 
# * snRNAseq: Binning performed across replicates within each sample, which 
#   always come from different regions (M, P, or A)
# 
# * RNAscope: Binning performed among replicates for each experiment (CIRCLE or 
#   STAR), followed by binning across experiments to collapse the disjoint cell
#   info across them.
#
# Saves the files:
#  * RNAscope data.frame: "df-rnascope-all_cell-prop-abund_dlpfc-ro1.rda"
#  * snRNAseq data.frame: "df-snrnaseq-all_cell-prop-abund_dlpfc-ro1.rda"
#

libv <- c("SingleCellExperiment", "SummarizedExperiment", "ggplot2")
sapply(libv, library, character.only = T)

#----------
# load data
#----------
# manage paths
base.path <- "."
save.dname <- "06_cellprop-integrated"
save.dpath <- file.path(base.path, "deconvo_method-paper", "outputs", save.dname)

# load sce
sce.fname <- "sef_mr-markers-k7-from-sce_dlpfc-ro1.rda"
sce.fpath <- file.path(base.path, "deconvo_method-paper", "outputs", 
                       "05_marker-gene-annotations", sce.fname)
sce <- get(load(sce.fpath))
# set vnames
donor.vname <- "BrNum"
ct.vname <- "cellType_broad_hc"

# load rnascope data
lcsv.fname <- "lcsv_halo.rda"
lcsv.fpath <- file.path(base.path, "deconvo_method-paper", 
                        "data", "halo", lcsv.fname)
if(!file.exists(lcsv.fpath)){
  fpath <- file.path(base.path, "Human_DLPFC_Deconvolution", "raw-data", "HALO")
  fnv <- list.files(fpath, pattern = ".*.csv$", recursive=T) # get long paths
  fnv <- fnv[!grepl("prelim", fnv)]; dirv <- c("CIRCLE", "STAR")
  lcsv <- lapply(dirv, function(diri){
    fnvi <- fnv[grepl(paste0(".*",diri,".*"), fnv)]
    lcsvi <- lapply(fnvi, function(fni){read.csv(file.path(fpath, fni))})
    names(lcsvi) <- gsub(".*\\/|\\.csv$", "", fnvi)
    return(lcsvi)
  })
  names(lcsv) <- dirv
  # save binary
  save(lcsv, file = lcsv.fpath)
} else{
  lcsv <- get(load(lcsv.fpath))
}


#-----------------
# helper functions
#-----------------
cellprop_from_halo <- function(csvi, samplab, region.str, expt = c("CIRCLE", "STAR"),
                               labexpt = list(CIRCLE = c("Endo", "Astro", "Inhib"),
                                              STAR = c("Excit", "Micro", "Oligo")),
                               labels = c("Endo" = "CLDN5", "Astro" = "GFAP",
                                          "Inhib" = "GAD1", "Excit" = "SLC17A7",
                                          "Micro" = "TMEM119", "Oligo" = "OLIG2")){
  # get labels for experiment
  labf <- labels[labexpt[[expt]]]
  # check cols
  cnvf <- which(grepl(paste(paste0("^", labf, "$"), collapse = "|"),
                      colnames(csvi)))
  dfi <- csvi[,cnvf]
  # relabel multi-mappers
  rsumv <- rowSums(dfi); which.mm <- which(rsumv > 1)
  dfi[which.mm,1] <- dfi[which.mm,2] <- dfi[which.mm,3] <- 0 
  # get proportions
  do.call(rbind, lapply(seq(length(labf)), function(ii){
    celli <- names(labf)[ii]; markeri <- labf[[ii]]
    num.celli <- length(which(dfi[,markeri]==1))
    data.frame(cell_type = celli, num_cells = num.celli, 
               prop_cells = num.celli/nrow(dfi), region = region.str, 
               sample_id = samplab, expt = expt)
  }))
}

# summarize tall df containing cell prop, counts etc.
get_dfmed_summary <- function(dfi, sampi){
  do.call(rbind, lapply(unique(dfi$cell_type), function(celli){
    dfii <- dfi[dfi$cell_type == celli,]
    data.frame(cell_type = celli,
               num_cells = paste(dfii$num_cells,collapse=";"),
               prop_cells = paste(dfii$prop_cells,collapse=";"),
               regions = paste(dfii$region, collapse = ";"),
               prop_median = median(dfii$prop_cells),
               num_median = median(dfii$num_cells),
               sample_id = sampi)
  }))
}

# get corr tests
get_lcor <- function(rni, sni, cnv = c("prop_median", "num_median")){
  lcor <- lapply(cnv, function(vari){cor.test(rni[,vari], sni[,vari])})
  names(lcor) <- cnv; return(lcor)
}

#-----------------------------
# get cell proportions, counts
#-----------------------------
# snrnaseq
sampv.sn <- unique(gsub("Br", "", sce[[donor.vname]]))
md <- colData(sce)
dfsn <- do.call(rbind, lapply(sampv.sn, function(sampi){
  filt.samp <- which(grepl(paste0("Br", sampi), md[,donor.vname]))
  mdi <- md[filt.samp,]
  dfi <- as.data.frame(table(mdi[, ct.vname]))
  dfi$prop <- dfi[,2]/sum(dfi[,2])
  colnames(dfi) <- c("cell_type", "num_cells", "prop_cells")
  dfi$sample_id <- sampi; dfi
}))
dfsn$assay <- "snRNAseq"

# rnascope slides
sampv <- unique(dfsn$sample_id)
dfrn <- do.call(rbind, lapply(c("CIRCLE", "STAR"), function(expt){
  lcsve <- lcsv[[expt]]
  do.call(rbind, lapply(sampv, function(sampi){
    message(sampi)
    str.patt <- paste0(".*_", sampi, "[A-Z].*")
    samp.filt <- grepl(str.patt, names(lcsve)) &  grepl("Final", names(lcsve))
    lcsvi <- lcsve[which(samp.filt)]
    if(length(lcsvi) > 0){
      do.call(rbind, lapply(seq(length(lcsvi)), function(ii){
        # get region string
        namei <- names(lcsvi)[ii]
        region.str <- gsub("[0-9]", "", strsplit(namei, "_")[[1]][3])
        region.str <- ifelse(region.str == "M", "middle",
                             ifelse(region.str == "A", "anterior", "posterior"))
        # make new cell prop df
        cellprop_from_halo(lcsvi[[ii]], sampi, region.str, expt)
      }))
    }
  }))
}))

# save datasets
# save.dpath <- "/users/smaden/"
rn.fname <- "df-rnascope-all_cell-prop-abund_dlpfc-ro1.rda"
sn.fname <- "df-snrnaseq-all_cell-prop-abund_dlpfc-ro1.rda"
save(dfrn, file = file.path(save.dpath, rn.fname))
save(dfsn, file = file.path(save.dpath, sn.fname))

#--------------------------
# cell proportion summaries
#--------------------------
# get summarized data
dffsn <- do.call(rbind, lapply(unique(dfsn$sample_id), function(sampi){
  get_dfmed_summary(dfsn[dfsn$sample_id==sampi,], sampi)
}))
dffrn <- do.call(rbind, lapply(unique(dfrn$sample_id), function(sampi){
  get_dfmed_summary(dfrn[dfrn$sample_id==sampi,], sampi)
}))

# filter overlapping samples, cell types
# relabel dfsn Micro.Oligo -> Micro
dfsn$cell_type <- as.character(dfsn$cell_type)
dfsn[dfsn$cell_type=="Micro.Oligo",]$cell_type <- "Micro"
sid.int <- intersect(dfsn$sample_id, dfrn$sample_id)
ct.int <- intersect(dfsn$cell_type, dfrn$cell_type)



which.int.sn <- dffsn$sample_id %in% sid.int &
  dffsn$cell_type %in% ct.int
which.int.rn <- dffrn$sample_id %in% sid.int &
  dffrn$cell_type %in% ct.int
dffsn <- dffsn[which.int.sn,]
dffrn <- dffrn[which.int.rn,]

#----------------
# total estimates
#----------------
# overall total est
ctv <- intersect(dffsn$cell_type, dffrn$cell_type)
dfm <- do.call(rbind, lapply(ctv, function(cti){
  rni <- dffrn[dffrn$cell_type==cti,]
  sni <- dffsn[dffsn$cell_type==cti,]
  int.sampv <- intersect(rni$sample_id, rni$sample_id)
  rni <- rni[rni$sample_id %in% int.sampv,]
  sni <- sni[sni$sample_id %in% int.sampv,]
  rni <- rni[order(match(rni$sample_id, sni$sample_id)),]
  cond <- identical(rni$sample_id, sni$sample_id)
  if(cond){
    data.frame(celltype = cti,
               sn.sum.num = sum(sni$num_median),
               sn.mean.num = mean(sni$num_median),
               sn.mean.prop = mean(sni$prop_median),
               rn.sum.num = sum(rni$num_median),
               rn.mean.num = mean(rni$num_median),
               rn.mean.prop = mean(rni$prop_median))
  }
}))
# corr across totals
corr.sum.num <- cor.test(dfm$sn.sum.num, dfm$rn.sum.num, method = "spearman")
corr.mean.prop <- cor.test(dfm$sn.mean.prop, dfm$rn.mean.prop, method = "spearman")
corr.mean.num <- cor.test(dfm$sn.mean.num, dfm$rn.mean.num, method = "spearman")

# get means across donors
ctv <- intersect(dffsn$cell_type, dffrn$cell_type)
dfmean <- do.call(rbind, lapply(ctv, function(cti){
  rni <- dffrn[dffrn$cell_type==cti,]
  sni <- dffsn[dffsn$cell_type==cti,]
  int.sampv <- intersect(rni$sample_id, rni$sample_id)
  rni <- rni[rni$sample_id %in% int.sampv,]
  sni <- sni[sni$sample_id %in% int.sampv,]
  rni <- rni[order(match(rni$sample_id, sni$sample_id)),]
  cond <- identical(rni$sample_id, sni$sample_id)
  if(cond){
    data.frame(celltype = cti,
               sn.num = mean(sni$num_median),
               sn.prop = mean(sni$prop_median),
               rn.num = mean(rni$num_median),
               rn.prop = mean(rni$prop_median))
  }
}))

#-----------------------------------------
# corr across cell types -- cell prop, num
#-----------------------------------------
ctv <- intersect(dffsn$cell_type, dffrn$cell_type)
df.cor <- do.call(rbind, lapply(ctv, function(cti){
  rni <- dffrn[dffrn$cell_type==cti,]
  sni <- dffsn[dffsn$cell_type==cti,]
  int.sampv <- intersect(rni$sample_id, rni$sample_id)
  rni <- rni[rni$sample_id %in% int.sampv,]
  sni <- sni[sni$sample_id %in% int.sampv,]
  rni <- rni[order(match(rni$sample_id, sni$sample_id)),]
  cond <- identical(rni$sample_id, sni$sample_id) &
    nrow(rni) > 2
  if(cond){
    corri.prop <- cor.test(rni$prop_median, sni$prop_median, method = "spearman")
    corri.num <- cor.test(rni$num_median, sni$num_median, method = "spearman")
    dfi <- matrix(c(corri.prop$estimate, corri.prop$p.value,
             corri.num$estimate, corri.num$p.value), nrow = 1)
    dfi <- as.data.frame(dfi)
  }
}))
colnames(df.cor) <- c("prop.rho", "prop.pval", "num.rho", "num.pval")
df.cor$prop.pbh <- p.adjust(df.cor$prop.pval, method = "BH")
df.cor$num.pbh <- p.adjust(df.cor$num.pval, method = "BH")

#---------------
# plot summaries
#---------------
# violin plots
ggplot(dfsn, aes(x = cell_type, y = prop_cells)) + geom_violin(draw_quantiles = 0.5)
ggplot(dfsn, aes(x = cell_type, y = num_cells)) + geom_violin(draw_quantiles = 0.5)








