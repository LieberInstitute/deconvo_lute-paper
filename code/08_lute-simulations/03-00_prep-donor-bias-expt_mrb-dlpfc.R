#!/usr/bin/env R

# Author: Sean Maden
#
# Prepare datasets for donor bias experiments using the multi-region brain 
# DLPFC dataset.
#
#

libv <- c("scuttle", "dplyr", "limma", "ggplot2", "ggforce", 
          "glmGamPoi", "sva", "DeconvoBuddies",
          "SingleCellExperiment", "SummarizedExperiment",
          "limma")
sapply(libv, library, character.only = TRUE)


#----------
# load data
#----------
# get save dpath
code.dname <- "08_lute-simulations"
proj.dname <- "deconvo_method-paper"
save.dpath <- file.path(proj.dname, "outputs", code.dname)
# get sce
sce.fname <- "sce-mrb_dlpfc.rda"
sce.fpath <- file.path(save.dpath, sce.fname)
sce <- get(load(sce.fpath))

#------------------
# format assay data
#------------------
for(ii in names(assays(sce))){
  assays(sce)[[ii]] <- as.matrix(assays(sce)[[ii]])}

# do limma removeBatchEffect
# mexpr <- assays(sce)[["counts"]]
# pheno <- data.frame(donor = sce[["donor"]], type = sce[["k2"]])
# mod <- model.matrix(~type, data = pheno)
# mi.adj <- removeBatchEffect(mexpr, batch = pheno$donor, covariates = mod)
# mi.adj[mi.adj < 0] <- 0 # convert negative values

# do combat-seq adjustment
# mexpr <- assays(sce)[["counts_ds"]]
# cnv <- colnames(sce)
# pheno <- data.frame(donor = sce[["donor"]], 
#                    type = sce[["k2"]])
# mod <- model.matrix(~type, data = pheno)
# mi.adj <- ComBat_seq(counts = mexpr, batch = pheno$donor, group = pheno$type)

# do combat adjustment
mexpr <- assays(sce)[["counts"]]
cnv <- colnames(sce)
pheno <- data.frame(donor = sce[["donor"]], 
                    type = sce[["k2"]])
mod <- model.matrix(~type, data = pheno)
mi.adj <- ComBat(dat = mexpr, batch = pheno$donor, mod = mod)
mi.adj[mi.adj < 0] <- 0 # convert negative values
# append to sce
assays(sce)[["counts_combat"]] <- mi.adj

# downsample on batches -- cell types
# note: downsample across donors within each time, then aggregate
assayname <- "counts"; typevar <- "k2"; batchvar <- "donor"
utypev <- unique(sce[[typevar]])

old.assay.name <- "counts_combat"
new.assay.name <- "counts_adj"
# iterate on types; downsample donors, then bind sce subsets/scef
sce <- do.call(cbind, lapply(utypev, function(ti){
  # filter sce
  sce.filt <- sce[[typevar]]==ti
  scef <- sce[,sce.filt]
  # get filtered data
  batchv <- scef[[batchvar]]
  mexpr <- assays(scef)[[old.assay.name]]
  # downsample
  mexpr.ds <- downsampleBatches(mexpr, batch = batchv)
  assays(scef)[[new.assay.name]] <- mexpr.ds 
  scef
}))

# append metadata
assays(sce) <- assays(sce)[names(assays(sce) %in% c("counts", "counts_adj"))]
metadata(sce)[["adj_method"]] <- "combat"
metadata(sce)[["ds_method"]] <- list(package = "scuttle",
                                     funct = "downsampleBatches",
                                     typevar = "k2", groupvar = "donor")

# get top markers
# get markers with mean ratios
celltypevar <- "k2"
sce <- logNormCounts(sce, assay.type = "counts_adj")
mr <- get_mean_ratio2(sce, assay_name = "logcounts", 
                      cellType_col = celltypevar)
# get top 20 markers from results
nmarker <- 20
typev <- unique(mr$cellType.target)
mr20 <- do.call(rbind, lapply(typev, function(typei){
  mr %>% filter(cellType.target == typei) %>% 
    arrange(rank_ratio) %>% top_n(n = nmarker)
}))
# store in sce
metadata(sce)[["k2.markers"]] <- list(all = mr, top20 = mr20)

#--------------------
# analyze dispersions
#--------------------
# plot dispersions
ldi <- sce_dispersion(sce, group.data = "k2", 
                      assayname = "counts_adj",
                      downsample = FALSE,
                      highlight.markers = mr20$gene,
                      hl.color = "red", point.alpha = 0.01,
                      hl.alpha = 1, ref.linecol = "black")
ldi$ggplot.dispersion
# append to sce data
metadata(sce)[["dispersion.counts_adj"]] <- ldi

# get dispersions by type
# compare dispersion coefficients from fitted neg binom 
# parse params
typevar <- "k2"; marker.name <- "k2_top20"
assay.name <- "counts_adj"
# get bg genes
num.genes.bg <- 1000
bg.name <- paste0("bg_", num.genes.bg)
genes.samplev <- sample(seq(nrow(sce)), num.genes.bg)
# get marker genes
genes.markerv <- metadata(sce)[["k2.markers"]][["top20"]]$gene
# define categories
catv <- c(unique(typev), "all") 
# get plot data
dfp <- do.call(rbind, lapply(catv, function(typei){
  # parse filter
  type.filt <- seq(ncol(sce))
  if(!typei == "all"){type.filt <- sce[[typevar]] == typei}
  scef <- sce[,type.filt]
  
  # get dispersions
  mexpr <- assays(scef)[[assay.name]]
  lglm.bg <- glm_gp(mexpr[genes.samplev,], on_disk = F)
  lglm.top20 <- glm_gp(mexpr[genes.markerv,], on_disk = F)
  
  # get plot data
  dfp1 <- data.frame(disp = lglm.bg$overdispersions)
  dfp1$marker.type <- bg.name
  dfp2 <- data.frame(disp = lglm.top20$overdispersions)
  dfp2$marker.type <- marker.name
  dfp <- rbind(dfp1, dfp2)
  dfp$celltype <- typei
  return(dfp)
}))
# set return list
ldisp <- list(dfp = dfp)

# boxplots at 3 zoom levels
ldisp[["ggbox"]] <- list()
ldisp[["ggbox"]][["zoom1"]] <- ggplot(dfp, aes(x = marker.type, y = disp)) + 
  geom_boxplot() + facet_wrap(~celltype)
ldisp[["ggbox"]][["zoom2"]] <- ggplot(dfp, aes(x = marker.type, y = disp)) + 
  geom_boxplot() + facet_wrap(~celltype) + ylim(0, 350)
ldisp[["ggbox"]][["zoom3"]] <- ggplot(dfp, aes(x = marker.type, y = disp)) + 
  geom_boxplot() + facet_wrap(~celltype) + ylim(0, 50)

# jitter plots at 3 zoom levels
ldisp[["ggjitter"]] <- list()
ldisp[["ggjitter"]][["zoom1"]] <- ggplot(dfp, aes(x = marker.type, y = disp)) + 
  geom_jitter(alpha = 0.5) + 
  stat_summary(geom = "crossbar", fun = "median", color = "red") + 
  facet_wrap(~celltype)
ldisp[["ggjitter"]][["zoom2"]] <- ggplot(dfp, aes(x = marker.type, y = disp)) + 
  geom_jitter(alpha = 0.5) + 
  stat_summary(geom = "crossbar", fun = "median", color = "red") + 
  facet_wrap(~celltype) + ylim(0, 350)
ldisp[["ggjitter"]][["zoom3"]] <- ggplot(dfp, aes(x = marker.type, y = disp)) + 
  geom_jitter(alpha = 0.5) + 
  stat_summary(geom = "crossbar", fun = "median", color = "red") + 
  facet_wrap(~celltype) + ylim(0, 50)

# append to sce
metadata(sce)[["dispersion.comparisons"]] <- ldisp

#---------------
# analyze anovas
#---------------
# ngene.sample = 5000
get_anova_df <- function(sce, ngene.sample = 2000, seed.num = 0,
                         assayv = c("counts", "counts_adj")){
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
    stat_summary(geom = "crossbar", fun = "median", color = "red") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_zoom(ylim = c(0, 10))
  # plot donor perc var
  lggj[["donor"]] <- ggplot(dfa.all, aes(x = assay, y = perc.var.donor)) +
    geom_jitter(alpha = 0.5) + ggtitle(title.str) +
    stat_summary(geom = "crossbar", fun = "median", color = "red") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_zoom(ylim = c(0, 2))
  # plot celltype.donor perc var
  lggj[["celltype.donor"]] <- 
    ggplot(dfa.all, aes(x = assay, y = perc.var.celltype.donor)) +
    geom_jitter(alpha = 0.5) + ggtitle(title.str) +
    stat_summary(geom = "crossbar", fun = "median", color = "red") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_zoom(ylim = c(0, 2))
  # plot residuals perc var  
  lggj[["residuals"]] <- 
    ggplot(dfa.all, aes(x = assay, y = perc.var.residuals)) +
    geom_jitter(alpha = 0.5) + ggtitle(title.str) +
    stat_summary(geom = "crossbar", fun = "median", color = "red") +
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

# append to sce object
metadata(sce)[["anova.analysis"]] <- lanova

#--------------
# save sce data
#--------------
sce.fname <- "sce_batch-adj_analysis-results_mrb-dlpfc.rda"
save(sce, file = file.path(save.dpath, sce.fname))

#-----------------------------
# get heatmaps and set objects
#-----------------------------
sce.fname <- "sce_batch-adj_analysis-results_mrb-dlpfc.rda"
sce <- get(load(file.path(save.dpath, sce.fname)))

ai <- "counts_adj"
mexpr <- assays(sce)[[ai]]
assays(sce)[[ai]] <- as.matrix(mexpr)
set1 <- set_from_sce(sce, assayname = ai, type.variable = "k2", get.group.stat = FALSE)




