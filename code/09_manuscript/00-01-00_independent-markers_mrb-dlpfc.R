#!/usr/bin/env R

#
#
#
#

libv <- c("lute", "scuttle", "dplyr", "limma", "ggplot2", "ggforce", 
          "glmGamPoi", "sva", "DeconvoBuddies",
          "SingleCellExperiment", "SummarizedExperiment",
          "limma")
sapply(libv, library, character.only = TRUE)


#----------
# load data
#----------
# get save dpath
code.dname <- "09_manuscript"
proj.dname <- "deconvo_method-paper"
save.dpath <- file.path(proj.dname, "outputs", code.dname)
# get sce
sce.fname <- "sce-mrb_dlpfc.rda"
sce.fpath <- file.path(save.dpath, sce.fname)
sce <- get(load(sce.fpath))

#-----------------------------------
# assign marker labels at variable k
#-----------------------------------
celltypevar <- "cellType"
table(sce[[celltypevar]])
# Astro    Excit_A    Excit_B    Excit_C    Excit_D    Excit_E    Excit_F    Inhib_A    Inhib_B    Inhib_C 
# 782        529        773        524        132        187        243        333        454        365 
# Inhib_D    Inhib_E    Inhib_F Macrophage      Micro      Mural      Oligo        OPC      Tcell 
# 413          7          8         10        388         18       5455        572          9

# define marker categories
sce[["k2"]] <- ifelse(grepl("^Excit.*|^Inhib.*", sce[[celltypevar]]), "neuron", "other")
sce[["k3"]] <- ifelse(grepl("^Excit.*", sce[[celltypevar]]), "Excit", 
                      ifelse(grepl("^Inhib.*", sce[[celltypevar]]), "Inhib", "other"))
sce[["k4"]] <- ifelse(grepl("^Excit.*", sce[[celltypevar]]), "Excit", 
                      ifelse(grepl("^Inhib.*", sce[[celltypevar]]), "Inhib", 
                             ifelse(grepl("^Oligo$", sce[[celltypevar]]), "Oligo", "other")))

#----------------------------------
# global params for plots, analyses
#----------------------------------
# coldata
batchvar <- "donor"
marker.typev <- c("k2", "k3", "k4") # for iterations

# assays for adjustments
assayname.unadj <- "counts"
assayname.marker <- "logcounts"
assayname.adj <- "counts_adj"

# number of top markers per type
nmarker <- 20

# dispersion mean-var plots
smooth.linecol <- "cornflowerblue"
hl.color <- "red"
ref.linecol <- "black"

# dispersion from fits
num.genes.bg <- 1000

#------------------------
# run adjustment workflow
#------------------------
# format assay
assays(sce)[[assayname.unadj]] <- as.matrix(assays(sce)[[assayname.unadj]])

# run the workflow
for(markeri in marker.typev){
  message("working on marker: ", markeri, "...")
  message("doing combat adj...")
  mexpr <- assays(sce)[[assayname.unadj]]
  cnv <- colnames(sce)
  pheno <- data.frame(donor = sce[["donor"]], type = sce[[markeri]])
  mod <- model.matrix(~type, data = pheno)
  mi.adj <- ComBat(dat = mexpr, batch = pheno$donor, mod = mod)
  message("converting negative values...")
  mi.adj[mi.adj < 0] <- 0 # convert negative values
  assays(sce)[[assayname.adj]] <- mi.adj
  
  message("performing downsampling...")
  # downsample -- by donor, within celltypes
  utypev <- unique(sce[[markeri]])
  sce <- do.call(cbind, lapply(utypev, function(ti){
    message("downsampling for type ", ti, "...")
    # filter sce
    sce.filt <- sce[[markeri]]==ti
    scef <- sce[,sce.filt]
    # get filtered data
    batchv <- scef[[batchvar]]
    mexpr <- assays(scef)[[assayname.adj]]
    # downsample
    mexpr.ds <- downsampleBatches(mexpr, batch = batchv)
    assays(scef)[[assayname.adj]] <- mexpr.ds 
    scef
  }))
  
  message("Getting top markers..")
  sce <- logNormCounts(sce, assay.type = assayname.adj)
  mr <- get_mean_ratio2(sce, assay_name = "logcounts", cellType_col = markeri)
  # get top N markers from results
  typev <- unique(mr$cellType.target)
  mrtop <- do.call(rbind, lapply(typev, function(typei){
    mr %>% filter(cellType.target == typei) %>% 
      arrange(rank_ratio) %>% top_n(n = nmarker)
  }))
  # store in sce
  metadata(sce)[["markers"]] <- list(all = mr, top = mrtop)
  
  # filter assays
  sce.filt <- names(assays(sce)) %in% c(assayname.unadj, assayname.adj)
  assays(sce) <- assays(sce)[sce.filt]
  
  # append metadata
  metadata(sce)[["adj_method"]] <- "combat"
  metadata(sce)[["ds_method"]] <- list(package = "scuttle",
                                       funct = "downsampleBatches",
                                       typevar = markeri, groupvar = "donor")
  
  # save
  fname <- paste0("sce_marker-adj-",markeri,"_mrb-dlpfc.rda")
  save(sce, file = file.path(save.dpath, fname))
  message("finished with marker ", markeri)
}

#------------------------
# run dispersion analyses
#------------------------
# run workflow
for(markeri in marker.typev){
  message("loading the data...")
  sce.fname <- paste0("sce_marker-adj-",markeri,"_mrb-dlpfc.rda")
  sce.fpath <- file.path(save.dpath, sce.fname)
  scei <- get(load(sce.fpath))
  
  message("doing mean-var analysis...")
  markerv <- metadata(sce)[["markers"]][["top"]]$gene
  # unadj
  disp.unadj <- sce_dispersion(scei, group.data = markeri, 
                               assayname = assayname.unadj,
                               highlight.markers = markerv,
                               hl.color = hl.color,
                               downsample = FALSE, point.alpha = 0.01,
                               ref.linecol = ref.linecol,
                               smooth.linecol = smooth.linecol)
  plot.fname <- paste0("ggpt-mean-var_",markeri,"-",assayname.unadj,"_ro1-dlpfc.jpg")
  jpeg(file.path(save.dpath, plot.fname), width = 8, height = 3.5, units = "in", res = 400)
  print(disp.unadj$ggplot.dispersion); dev.off()
  # adj
  disp.adj <- sce_dispersion(scei, group.data = markeri, 
                             assayname = assayname.adj,
                             highlight.markers = mrtop$gene, 
                             hl.color = hl.color,
                             downsample = FALSE, point.alpha = 0.01,
                             ref.linecol = ref.linecol,
                             smooth.linecol = smooth.linecol)
  plot.fname <- paste0("ggpt-mean-var_",markeri,"-",assayname.adj,"_ro1-dlpfc.jpg")
  jpeg(file.path(save.dpath, plot.fname), width = 8, height = 3.5, units = "in", res = 400)
  print(disp.adj$ggplot.dispersion); dev.off()
  # append results to metadata
  metadata(scei)[["disp_mean-var_unadj"]] <- disp.unadj
  metadata(scei)[["disp_mean-var_adj"]] <- disp.adj
  
  message("analyzing dispersions from nb fits...") 
  # parse params
  typevar <- "k2"; marker.name <- "k2_top20"
  assay.name <- "counts_adj"
  # get bg genes
  
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
  ldisp[["ggbox"]] <- ggplot(dfp, aes(x = marker.type, y = disp)) + 
    geom_boxplot() + facet_wrap(~celltype)
  # jitter plots at 3 zoom levels
  ldisp[["ggjitter"]] <- ggplot(dfp, aes(x = marker.type, y = disp)) + 
    geom_jitter(alpha = 0.5) + 
    stat_summary(geom = "crossbar", fun = "median", color = "red") + 
    facet_wrap(~celltype)
  # append to sce
  metadata(sce)[["dispersion_fits"]] <- ldisp
  
  # save
  save(sce, file = file.path(save.dpath, sce.fname))
  
}