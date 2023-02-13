#!/usr/bin/env R

#
#
#
#

libv <- c("lute", "scuttle", "dplyr", "limma", "ggplot2", "ggforce", "gridExtra",
          "glmGamPoi", "sva", "DeconvoBuddies", "SingleCellExperiment", "limma",
          "SummarizedExperiment")
sapply(libv, library, character.only = TRUE)

#---------------
# manage paths
#---------------
# get save dpath
code.dname <- "09_manuscript"
proj.dname <- "deconvo_method-paper"
save.dpath <- file.path(proj.dname, "outputs", code.dname)

#----------------------------------
# global params for plots, analyses
#----------------------------------
# stem for new filenames
proj.stem <- "mrb-dlpfc"

# sce filename
sce.fname <- "sce-mrb_dlpfc.rda"

# coldata
celltypevar <- "cellType"
batchvar <- "donor"
marker.typev <- c("k2", "k3", "k4") # for iterations

# assays for adjustments
assayname.unadj <- "counts"
assayname.marker <- "logcounts"
assayname.adj <- "counts_adj"

# number of top markers per type
nmarker <- 20

# scale filter
scale.thresh <- -0.9

# dispersion mean-var plots
smooth.linecol <- "cornflowerblue"
hl.color <- "red"
ref.linecol <- "black"

# dispersion from fits
num.genes.bg <- 1000

#----------
# load data
#----------
# get sce
sce.fpath <- file.path(save.dpath, sce.fname)
sce <- get(load(sce.fpath))

#------------------
# filter cell types
#------------------
sce.filt <- sce[[celltypevar]] %in% c("Macrophage", "Mural", "Tcell")
sce <- sce[,!sce.filt]

table(sce[[celltypevar]])
# Astro    Excit_A    Excit_B    Excit_C    Excit_D    Excit_E    Excit_F    Inhib_A    Inhib_B 
# 782        529        773        524        132        187        243        333        454 
# Inhib_C    Inhib_D    Inhib_E    Inhib_F Macrophage      Micro      Mural      Oligo        OPC 
# 365        413          7          8          0        388          0       5455        572 
# Tcell 
# 0 

#-----------------------------------
# assign marker labels at variable k
#-----------------------------------
# define marker categories
sce[["k2"]] <- ifelse(grepl("^Excit.*|^Inhib.*", sce[[celltypevar]]), 
                      "neuron", "glial")
sce[["k3"]] <- ifelse(grepl("^Excit.*", sce[[celltypevar]]), "Excit", 
                      ifelse(grepl("^Inhib.*", sce[[celltypevar]]), 
                             "Inhib", "glial"))
sce[["k4"]] <- ifelse(grepl("^Excit.*", sce[[celltypevar]]), "Excit", 
                      ifelse(grepl("^Inhib.*", sce[[celltypevar]]), "Inhib", 
                             ifelse(grepl("^Oligo$", sce[[celltypevar]]), "Oligo", 
                                    "non_oligo_glial")))

#--------------
# format assays
#--------------
assays(sce)[[assayname.unadj]] <- as.matrix(assays(sce)[[assayname.unadj]])

#------------------------
# run adjustment workflow
#------------------------
# run the workflow
for(markeri in marker.typev){
  message("working on marker: ", markeri, "...")
  message("filtering low-abundance donors by type...")
  cd <- colData(sce)
  dfm <- as.data.frame(table(cd[,markeri]))
  type.filt <- as.character(dfm[dfm[,2]==min(dfm[,2]),1])
  cdf <- cd[cd[,markeri] == type.filt,]
  dft <- as.data.frame(table(cdf[,batchvar], cdf[,markeri]))
  scalev <- scale(dft$Freq)[,1]
  which.filt <- which(scalev <= scale.thresh)
  sample.filt <- as.character(dft[which.filt,1])
  cells.remove <- rownames(cd[cd[,batchvar] %in% sample.filt,])
  sce <- sce[,!colnames(sce) %in% cells.remove]
  md.rm <- list(type.filt = type.filt,
                scale.thresh = scale.thresh,
                num.cells = length(cells.remove),
                batch.id = sample.filt)
  message("removed ", length(cells.remove),
          " cells from ", length(sample.filt),
          " batches for cell type: '", type.filt, "'.")
  
  message("doing combat adj...")
  mexpr <- assays(sce)[[assayname.unadj]]
  cnv <- colnames(sce)
  pheno <- data.frame(donor = sce[[batchvar]], type = sce[[markeri]])
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
  metadata(sce)[["scale_filter"]] <- md.rm
  metadata(sce)[["adj_method"]] <- "combat"
  metadata(sce)[["ds_method"]] <- list(package = "scuttle",
                                       funct = "downsampleBatches",
                                       typevar = markeri, groupvar = "donor")
  
  # save
  fname <- paste0("sce_marker-adj-",markeri,"_",proj.stem,".rda")
  save(sce, file = file.path(save.dpath, fname))
  message("finished with marker ", markeri)
}

#----------------------
# get marker expression
#----------------------
marker.typev <- c("k2", "k3", "k4")

lscef <- lapply(marker.typev, function(markeri){
  message("loading the data...")
  sce.fname <- paste0("sce_marker-adj-",markeri,"_",
                      proj.stem,".rda")
  sce.fpath <- file.path(save.dpath, sce.fname)
  scei <- get(load(sce.fpath))
  
  # get marker expr
  mr <- metadata(sce)[["markers"]][["top"]] # top markers
  scef <- scei[mr$gene,]
  return(scef)
})
names(lscef) <- marker.typev

# save
fname <- paste0("list-scef_markers-k2-k3-k4_",proj.stem,".rda")
save(lscef, file = file.path(save.dpath, fname))

#------------------------
# run dispersion analyses
#------------------------
set.seed(0) # for sampling from bg

# run workflow
for(markeri in marker.typev){
  message("loading the data...")
  sce.fname <- paste0("sce_marker-adj-",markeri,"_",proj.stem,".rda")
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
  plot.fname <- paste0("ggpt-mean-var_",markeri,"-",assayname.unadj,"_",proj.stem,".jpg")
  jpeg(file.path(save.dpath, plot.fname), width = 8, height = 3.5, units = "in", res = 400)
  print(disp.unadj$ggplot.dispersion); dev.off()
  # adj
  disp.adj <- sce_dispersion(scei, group.data = markeri, 
                             assayname = assayname.adj,
                             highlight.markers = markerv, 
                             hl.color = hl.color,
                             downsample = FALSE, point.alpha = 0.01,
                             ref.linecol = ref.linecol,
                             smooth.linecol = smooth.linecol)
  plot.fname <- paste0("ggpt-mean-var_",markeri,"-",
                       gsub("_", "-", assayname.adj),
                       "_ro1-dlpfc.jpg")
  jpeg(file.path(save.dpath, plot.fname), width = 8, height = 3.5, units = "in", res = 400)
  print(disp.adj$ggplot.dispersion); dev.off()
  # append results to metadata
  metadata(scei)[["disp_mean-var_unadj"]] <- disp.unadj
  metadata(scei)[["disp_mean-var_adj"]] <- disp.adj
  
  message("analyzing dispersions from nb fits...") 
  # labels
  marker.name <- "top_markers"
  bg.name <- paste0("bg_", num.genes.bg)
  # get bg genes
  genes.samplev <- sample(seq(nrow(sce)), num.genes.bg)
  # get marker genes
  genes.markerv <- metadata(sce)[["markers"]][["top"]]$gene
  # define categories
  typev <- sce[[markeri]]
  catv <- c(unique(typev), "all") 
  # get plot data
  assayv <- c(assayname.unadj, assayname.adj)
  dfp <- do.call(rbind, lapply(catv, function(typei){
    # parse filter
    type.filt <- seq(ncol(scei))
    if(!typei == "all"){type.filt <- scei[[markeri]] == typei}
    scef <- scei[,type.filt]
    # iterate on assays
    dfpi <- do.call(rbind, lapply(assayv, function(assayi){
      # get dispersions
      mexpr <- assays(scef)[[assayi]]
      lglm.bg <- glm_gp(mexpr[genes.samplev,], on_disk = F)
      lglm.top <- glm_gp(mexpr[genes.markerv,], on_disk = F)
      # get plot data
      dfp1 <- data.frame(disp = lglm.bg$overdispersions)
      dfp1$marker.type <- bg.name
      dfp1$marker <- rownames(mexpr[genes.samplev,])
      dfp2 <- data.frame(disp = lglm.top$overdispersions)
      dfp2$marker.type <- marker.name
      dfp2$marker <- rownames(mexpr[genes.markerv,])
      dfp <- rbind(dfp1, dfp2)
      dfp$assay <- assayi
      return(dfp)
    }))
    dfpi$celltype <- typei
    return(dfpi)
  }))
  # set return list
  ldisp <- list(dfp = dfp)
  # append to sce
  metadata(sce)[["dispersion_fits"]] <- ldisp
  
  message("saving updated sce object...")
  save(sce, file = file.path(save.dpath, sce.fname))
  message("finished with marker group ", markeri, ".")
}

# make new dispersion boxplots
for(markeri in marker.typev){
  message("working on marker ", markeri, "...")
  message("loading the data...")
  sce.fname <- paste0("sce_marker-adj-",markeri,"_",proj.stem,".rda")
  sce.fpath <- file.path(save.dpath, sce.fname)
  scei <- get(load(sce.fpath))
  
  # get plot data
  dfp <- metadata(scei)[["dispersion_fits"]][["dfp"]]
  
  # boxplots with manual ylim, facet on celltype;assay
  plot1 <- ggplot(dfp, aes(x = marker.type, y = disp)) + 
    theme_bw() + geom_boxplot() + ylim(0, 100) +
    facet_wrap(~celltype+assay, nrow = 2) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  fname <- paste0("ggbox_bg1k-vs-topmarkers-",markeri,"_",proj.stem,".jpg")
  jpeg(file.path(save.dpath, fname), width = 8, height = 4, 
       units = "in", res = 400)
  print(plot1); dev.off()
  
  # boxplots with zooms, use facet zoom
  for(typei in unique(dfp$celltype)){
    ymax.zoom <- 15
    yaxis.title <- paste0(paste0(rep(" ", 15), collapse = ""), "Dispersion")
    fname <- paste0("ggbox-",typei,
                    "_bg1k-vs-topmarkers-",markeri,
                    "_",proj.stem,".jpg")
    title.str <- paste0("Unadj. counts, ", typei)
    
    dfp.filt <- dfp$assay=="counts" & dfp$celltype==typei
    dfpf <- dfp[dfp.filt,]
    plot2.1 <- ggplot(dfpf, aes(x = marker.type, y = disp)) + 
      theme_bw() + geom_boxplot() + facet_zoom(ylim = c(0, ymax.zoom)) +
      ggtitle(title.str) +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_blank())
    title.str <- paste0("Adj. counts, ", typei)
    
    dfp.filt <- dfp$assay=="counts_adj" & dfp$celltype==typei
    dfpf <- dfp[dfp.filt,]
    plot2.2 <- ggplot(dfpf, aes(x = marker.type, y = disp)) + 
      theme_bw() + geom_boxplot() + facet_zoom(ylim = c(0, ymax.zoom)) +
      ggtitle(title.str) + 
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1))
    
    jpeg(file.path(save.dpath, fname), width = 5, height = 3, 
         units = "in", res = 400)
    grid.arrange(plot2.1, plot2.2, nrow = 2, 
                 bottom = "Genes",
                 layout_matrix = matrix(c(1,1,2,2,2), ncol = 1),  
                 left = yaxis.title)
    dev.off()
  }
  
  # scatter plots of unadj vs. adj marker disp
  ymax <- 20; xmax <- 20; ptcol <- "red"
  for(typei in unique(dfp$celltype)){
    message("working on scatterplot for type ", typei, "...")
    dfpi <- dfp[dfp$celltype==typei,]
    dfpi <- dfpi[dfpi$marker.type=="top_markers",]
    dfp1 <- dfpi[dfpi$assay=="counts",]
    dfp2 <- dfpi[dfpi$assay=="counts_adj",]
    dfp.new <- data.frame(unadj = dfp1$disp, adj = dfp2$disp, 
                          marker = dfp1$marker)
    title.str <- paste0(typei, " expression at top markers (", 
                        length(unique(dfp.new$marker)), " genes)")
    ggpt <- ggplot(dfp.new, aes(x = unadj, y = adj)) + theme_bw() +
      geom_point(alpha = 0.5, color = ptcol) + 
      geom_abline(slope = 1, intercept = 0) + 
      facet_zoom(ylim = c(0, ymax), xlim = c(0, xmax)) +
      xlab("Unadjusted dispersion") + ylab("Adjusted dispersion") + 
      ggtitle(title.str)
    
    fname <- paste0("ggpt-markers_", markeri, "-", typei,"_",proj.stem,".jpg")
    jpeg(file.path(save.dpath, fname), 
         width = 5, height = 3, units = "in", res = 400)
    print(ggpt); dev.off()
  }
  
  message("finished with marker ", markeri, ".")
}
