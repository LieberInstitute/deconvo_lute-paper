#!/usr/bin/env R

# Author: Sean Maden
#
# Prepare datasets for donor bias experiments using the multi-region brain 
# DLPFC dataset.
#
#

libv <- c("scuttle", "dplyr", "ggplot2", "ggforce", 
          "glmGamPoi", "sva", "DeconvoBuddies",
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
# dfm.fname <- "markers-k2_db-mr2_sce-dlpfc-mrb.rda"
# dfm <- get(load(file.path(save.dpath, dfm.fname)))
# get sce
sce.fname <- "sce-mrb_dlpfc.rda"
sce.fpath <- file.path(save.dpath, sce.fname)
sce <- get(load(sce.fpath))

#------------------
# format assay data
#------------------
for(ii in names(assays(sce))){
  assays(sce)[[ii]] <- as.matrix(assays(sce)[[ii]])}

#---------------------
# perform downsampling
#---------------------
mct <- assays(sce)[["counts"]]
# mlc <- assays(sce)[["logcounts"]]

# downsample on proportions vector
# assays(sce)[["counts_ds-mat"]] <- mexpr_downsample(mct)
# assays(sce)[["logcounts_ds-mat"]] <- mexpr_downsample(mlc)
# downsample on batches -- donors
# assays(sce)[["counts_ds-donor"]] <- downsampleBatches(mct, batch = sce[["donor"]])
# assays(sce)[["logcounts_ds-donor"]] <- downsampleBatches(mlc, batch = sce[["donor"]])

# downsample on batches -- cell types
# note: downsample across donors within each time, then aggregate
assayname <- "counts"; typevar <- "k2"; batchvar <- "donor"
utypev <- unique(sce[[typevar]])
sce <- do.call(cbind, lapply(utypev, function(ti){
  sce.filt <- sce[[typevar]]==ti; scef <- sce[,sce.filt]
  batchv <- scef[[batchvar]]; mexpr <- assays(scef)[[assayname]]
  mexpr.ds <- downsampleBatches(mexpr, batch = batchv)
  assays(scef)[[paste0(assayname, "_ds")]] <- mexpr.ds 
  scef
}))

# save
sce.fname <- "sce_ds-assays_mrb-dlpfc.rda"
save(sce, file = file.path(save.dpath, sce.fname))

#---------------
# perform combat
#---------------
assayv <- c("counts", "counts_ds")
for(ai in assayv){
  message("Working on assay: ", ai)
  mi <- assays(sce)[[ai]]; cnv <- colnames(sce)
  pheno <- data.frame(donor = sce[["donor"]], type = sce[["k2"]])
  mod <- model.matrix(~type, data = pheno)
  batch <- pheno$donor
  mi.adj <- ComBat(dat = mi, batch = batch, mod = mod, 
                   par.prior = TRUE, prior.plots = FALSE)
  assays(sce)[[paste0(ai, "_combat")]] <- mi.adj
}

# convert negative values
assayname <- "counts_ds_combat"
me <- assays(sce)[[assayname]]
me[me < 0] <- 0
assays(sce)[[assayname]] <- me

# filter sce -- cleanup unused assays
#assays.remove <- paste0(rep(c("counts", "logcounts"), 2), 
#                        rep(c("_ds-mat", "_ds-donor"), each = 2))
assays.remove <- c("logcounts", "counts_ds")
assays(sce) <- assays(sce)[!names(assays(sce)) %in% assays.remove]

# save
sce.fname <- "sce_combat-ds-assays_mrb-dlpfc.rda"
save(sce, file = file.path(save.dpath, sce.fname))

#-------------------
# get top 20 markers
#-------------------
# get markers with mean ratios
sce <- logNormCounts(sce, assay.type = "counts_ds_combat")
mr <- get_mean_ratio2(sce, assay_name = "logcounts", cellType_col = "k2")

# get top 20 markers from results
nmarker <- 20
typev <- unique(mr$cellType.target)
mr20 <- do.call(rbind, lapply(typev, function(typei){
  mr %>% filter(cellType.target == typei) %>% 
    arrange(rank_ratio) %>% top_n(n = nmarker)
}))

# store in sce
metadata(sce)[["k2.markers"]] <- list(all = mr, top20 = mr20)

# save sce with marker data
sce.fname <- "sce_combat-ds-assays_mrb-dlpfc.rda"
save(sce, file = file.path(save.dpath, sce.fname))

#-------------------
# analysis functions
#------------------

library(limma)

# ngene.sample = 5000
get_anova_df <- function(sce, ngene.sample = NULL, seed.num = 0,
                         assayv = c("counts", "counts_ds_combat")){
  set.seed(seed.num)
  lr <- lgg <- lggj <- lggpt <- list()
  if(is(ngene.sample, "NULL")){
    sampv <- seq(nrow(sce))
  } else{
    sampv <- sample(seq(nrow(sce)), ngene.sample, replace = FALSE)  
  }
  dfa.all <- do.call(rbind, lapply(assayv, function(ai){
    message("working on assay: ", ai, "...")
    mi <- assays(sce)[[ai]];
    mi <- mi[sampv,] # get random subset
    maxv <- rowMaxs(mi); mi <- mi[maxv > 0,] # filter all-zeros
    num.na <- apply(mi, 1, function(ri){length(which(is.na(ri)))}) 
    mi <- mi[which(num.na == 0),] # filter nas
    dfi <- data.frame(expr = mi[1,], celltype = sce[["k2"]],
                      donor = sce[["donor"]])
    dfa.mi <- do.call(rbind, lapply(seq(nrow(mi)), function(ii){
      message(ii)
      dfi$expr <- mi[ii,]
      
      avi <- aov(formula = expr ~ celltype * donor, data = dfi)
      
      # old method:
      # lmi <- lm(expr ~ celltype * donor, data = dfi)
      # avi <- anova(lmi)
      
      # use limma
      # dmat <- model.matrix(~ celltype * donor, data = dfi)
      # avi <- limma::lmFit(object = mi, design = dmat)
      
      # perc.var <- 100*avi$`Sum Sq`/sum(avi$`Sum Sq`)
      ssqv <- summary(avi)[[1]][[2]]
      perc.var <- 100*ssqv/sum(ssqv)
      data.frame(perc.var.celltype = perc.var[1],
                 perc.var.donor = perc.var[2],
                 perc.var.celltype.donor = perc.var[3],
                 perc.var.residuals = perc.var[4],
                 marker = rownames(mi)[ii])
    }))
    message("finished anovas")
    dfa.mi <- dfa.mi[!is.na(dfa.mi[,1]),] # filter nas
    dfa.mi$assay <- ai; dfa.mi
  }))
  # plots
  num.genes <- length(unique(dfa.all$marker))
  title.str <- paste0("Perc. var. (",num.genes," genes)")
  # plot celltype perc var
  lggj[["celltype"]] <- 
    ggplot(dfa.all, aes(x = assay, y = perc.var.celltype)) +
    geom_jitter(alpha = 0.5) + ggtitle(title.str) +
    stat_summary(geom = "crossbar", fun = "mean", color = "red") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_zoom(ylim = c(0, 10))
  # plot donor perc var
  lggj[["donor"]] <- ggplot(dfa.all, aes(x = assay, y = perc.var.donor)) +
    geom_jitter(alpha = 0.5) + ggtitle(title.str) +
    stat_summary(geom = "crossbar", fun = "mean", color = "red") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_zoom(ylim = c(0, 2))
  # plot celltype.donor perc var
  lggj[["celltype.donor"]] <- 
    ggplot(dfa.all, aes(x = assay, y = perc.var.celltype.donor)) +
    geom_jitter(alpha = 0.5) + ggtitle(title.str) +
    stat_summary(geom = "crossbar", fun = "mean", color = "red") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_zoom(ylim = c(0, 2))
  # plot residuals perc var  
  lggj[["residuals"]] <- 
    ggplot(dfa.all, aes(x = assay, y = perc.var.residuals)) +
    geom_jitter(alpha = 0.5) + ggtitle(title.str) +
    stat_summary(geom = "crossbar", fun = "mean", color = "red") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_zoom(ylim = c(90, 100))
  
  # scatter plot of perc var, within assay
  get_ggpt <- function(a1 = "counts", a2 = "counts_ds_combat", 
                       regex.cnv = "^perc\\.var\\..*"){
    # get dfs
    df1 <- dfa.all[dfa.all$assay==a1,]
    df2 <- dfa.all[dfa.all$assay==a2,]
    # match dfs
    df1 <- df1[!duplicated(df1$marker),]
    df2 <- df2[!duplicated(df2$marker),]
    markerv.int <- intersect(df1$marker, df2$marker)
    df1 <- df1[df1$marker %in% markerv.int,]
    df2 <- df2[df2$marker %in% markerv.int,]
    df2 <- df2[order(match(df2$marker, df1$marker)),]
    # get plots if match successful
    cond <- identical(as.character(df2$marker), as.character(df1$marker))
    if(cond){
      # get colnames to plot
      cnv <- colnames(df1); cnv <- cnv[grepl(regex.cnv, cnv)]
      # plot each colname
      lggi <- lapply(cnv, function(ci){
        dfp <- data.frame(pv1 = df1[,ci], pv2 = df2[,ci])
        title.str <- paste0("Percent var. ", gsub(".*\\." , "", ci), 
                            "\n(",nrow(dfp)," genes)")
        ggplot(dfp, aes(x = pv1, y = pv2)) + geom_point(alpha = 0.5) + 
          geom_abline(slope = 1, intercept = 0, col = "black") +
          ggtitle(title.str) + xlab(a1) + ylab(a2)
      })
      names(lggi) <- cnv
      return(lggi)
    } else{
      stop("Error, couldn't match markers across provided assays.")
    }
    return(NULL)
  }
  lggpt <- get_ggpt()
  lgg <- list(gg.jitter = lggj, gg.pt = lggpt)
  lr <- list(dfa = dfa.all, lgg = lgg)
  return(lr)
}

# get anova analysis
lanova <- get_anova_df(sce)

# save
fname <- paste0("anova-results_mrb-dlpfc.rda")
save(lanova, file = file.path(save.dpath, fname))

#----------------------------
# get dispersion test results
#----------------------------
mr.hl <- metadata(sce)[["k2.markers"]][["top20"]]
markerv <- unique(mr.hl$gene)
assayv <- c("counts", "counts_ds_combat")

lld <- lapply(assayv, function(ai){
  message("working on assay: ", ai, "...")
  ldi <- sce_dispersion(sce, group.data = "k2", assayname = ai,
                        highlight.markers = mr.hl$gene, downsample = FALSE,
                        hl.color = "red", point.alpha = 0.01,
                        hl.alpha = 1, ref.linecol = "black")
  lf.all <- mexpr_nbcoef(as.matrix(assays(scef)[[ai]]))
  lf.k2 <- mexpr_nbcoef(as.matrix(assays(sce[mr.hl$gene,])[[ai]]))
  #ldi[["nb.fit"]] <- list(all.sub = lf.all, markers.k2 = lf.k2)
  return(ldi)
})
names(lld) <- assayv

# save
fname <- "dispersion-results_mrb-dlpfc.rda"
save(lld, file = file.path(save.dpath, fname))

#-------------------------
# get overdispersion plots
#-------------------------
# counts
mct <- assays(sce)[]
lf.all <- mexpr_nbcoef(mf1)




