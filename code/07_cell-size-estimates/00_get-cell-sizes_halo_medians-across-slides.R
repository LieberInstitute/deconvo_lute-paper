#!/usr/bin/env R

# Get the cell sizes from HALO outputs.
#
# Reads in HALO csvs with data.table::fread, using check.names = F to avoid 
# adding "." to column names.
#

libv <- c("ggplot2", "gridExtra", "ggpubr", "data.table")
sapply(libv, library, character.only = T)

#----------
# load data
#----------
# manage paths
base.path <- "."
save.dname <- "07_cell-size-estimates"
subdir.name <- "Algorithm_Check_20220920"
out.fnstem <- "algocheck-20220920"

save.dpath <- file.path(base.path, "deconvo_method-paper", "outputs", save.dname)

# load rnascope data
lcsv.fname <- "lcsv_halo.rda"
lcsv.fpath <- file.path(base.path, "deconvo_method-paper", "data", "halo", lcsv.fname)

if(!file.exists(lcsv.fpath)){
  fpath <- file.path(base.path, "Human_DLPFC_Deconvolution", "raw-data", "HALO",
                     subdir.name)
  fnv <- list.files(fpath, pattern = ".*.csv$", recursive=T) # get long paths
  fnv <- fnv[!grepl("prelim", fnv)]; dirv <- c("CIRCLE", "STAR")
  lcsv <- lapply(dirv, function(diri){
    fnvi <- fnv[grepl(paste0(".*",diri,".*"), fnv)]
    lcsvi <- lapply(fnvi, function(fni){
        data.table::fread(file.path(fpath, fni), sep = ",", check.names = F)
      }
    )
    names(lcsvi) <- gsub(".*\\/|\\.csv$", "", fnvi)
    return(lcsvi)
  })
  names(lcsv) <- dirv
  # save binary
  save(lcsv, file = lcsv.fpath)
} else{
  lcsv <- get(load(lcsv.fpath))
}

#-----------
# set params
#-----------
labels = c("Endo" = "CLDN5", "Astro" = "GFAP",
           "Inhib" = "GAD1", "Excit" = "SLC17A7",
           "Micro" = "TMEM119", "Oligo" = "OLIG2")

cnv.size <- c("nucleus_area" = "Nucleus Area (µm²)",
              "nucleus_perimenter" = "Nucleus Perimeter (µm)", 
              "akt3_expr" = "AKT3 (Opal 570) Copies")

# get flattened, filtered list
lc1 <- lapply(lcsv[[1]], function(ii){ii})
lc2 <- lapply(lcsv[[2]], function(ii){ii})
lc <- c(lc1, lc2)

#----------------------
# get all sizes by cell
#----------------------
dfs.cell <- do.call(rbind, lapply(seq(length(lc)), function(ii){
  message(ii)
  csvi <- as.data.frame(lc[[ii]])
  csv.fname <- names(lc)[ii]
  cnv <- colnames(csvi)
  # get cell type label
  cnv.type <- colnames(csvi)[grepl(".*Positive$", colnames(csvi))]
  cnv.type <- cnv.type[grepl(paste0(labels, collapse = "|"), cnv.type)]
  do.call(rbind, lapply(cnv.type, function(typei){
    message(typei)
    # get filter terms
    which.type <- which(csvi[,typei]==1) # celltype filter
    csvif <- csvi[which.type,]
    # get size variable
    csizev <- as.character(cnv.size)
    cnv.in <- intersect(colnames(csvif), csizev)
    csvif <- csvif[,cnv.in]
    for(ci in csizev[!csizev %in% colnames(csvif)]){
      csvif[,ci] <- NA
    }
    csvif <- csvif[,csizev]
    csvif[,"type"] <- names(labels[labels==gsub(" .*", "", typei)])
    csvif[,"csv.fname"] <- csv.fname; 
    csvif
  }))
}))

for(c in seq(3)){dfs.cell[,c] <- as.numeric(dfs.cell[,c])}
dfs.cell$slide <- unlist(lapply(dfs.cell$csv.fname, function(fni){
  unlist(strsplit(fni, "_"))[2]
}))
dfs.cell$position <- gsub("[0-9]", "", dfs.cell$slide)
dfs.cell$donor <- gsub("[A-Z]", "", dfs.cell$slide)

# save 
save.fname <- "df-csize-cells_halo.rda"
save(dfs.cell, file = file.path(save.dpath, save.fname))

# get plot data
# 
dfsc <- get(load(file.path(save.dpath, save.fname)))
cnv <- colnames(dfsc)
datv <- c(dfsc[,1], dfsc[,2], dfsc[,3])
dfp <- data.frame(value = datv, type = rep(dfsc$type, 3),
                  method = c(rep("nuc_area", nrow(dfsc)),
                             rep("nuc_perim", nrow(dfsc)),
                             rep("akt3_expr", nrow(dfsc))),
                  slide = c(rep(dfsc$slide, 3)))

# violin plots
ggvp <- ggplot(dfp, aes(x = type, y = value, color = type)) +
  geom_violin(draw_quantiles = 0.5) + theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())
ggvp <- ggvp + facet_wrap(~method)
# save new vp
plot.fname <- paste0("ggviolin-cells_csize-3vars_halo_",
                     out.fnstem,".jpg")
jpeg(file.path(save.dpath, plot.fname), 
     width = 5, height = 1.8, units = "in", res = 400)
ggvp; dev.off()
# new vp with jitter
plot.fname <- paste0("ggviolin-withjitter-cells_csize-3vars_halo_",
                     out.fnstem,".jpg")
jpeg(file.path(save.dpath, plot.fname), 
     width = 5, height = 1.8, units = "in", res = 400)
ggvp + geom_jitter(alpha = 0.5, size = 0.1); dev.off()

# plots by slide
# boxplots of values by slide
ggbp <- ggplot(dfp, aes(x = slide, y = value, color = type)) + geom_boxplot() +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggbp <- ggbp + facet_wrap(~method*type, nrow = 3)
# new boxplots
plot.fname <- paste0("ggboxplot-cells_csize-3vars_halo_",
                     out.fnstem,".jpg")
jpeg(file.path(save.dpath, plot.fname), 
     width = 9, height = 6.5, units = "in", res = 400)
ggbp; dev.off()
# new violin plots
ggvp <- ggplot(dfp, aes(x = slide, y = value, color = type)) + 
  geom_violin(draw_quantiles = 0.5) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggvp <- ggvp + facet_wrap(~method*type, nrow = 3)
plot.fname <- paste0("ggviolin-cells_csize-3vars_halo_",
                     out.fnstem,".jpg")
jpeg(file.path(save.dpath, plot.fname), 
     width = 9, height = 6.5, units = "in", res = 400)
ggvp; dev.off()
# new violin plots -- with jitter
ggvp <- ggplot(dfp, aes(x = slide, y = value, color = type)) + 
  geom_violin(draw_quantiles = 0.5) + geom_jitter(alpha = 0.1) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggvp <- ggvp + facet_wrap(~method*type, nrow = 3)
plot.fname <- paste0("ggviolin-jitter-cells_csize-3vars_halo_",
                     out.fnstem,".jpg")
jpeg(file.path(save.dpath, plot.fname), 
     width = 9, height = 6.5, units = "in", res = 400)
ggvp; dev.off()

# medians by slide
dfp.med <- do.call(rbind, lapply(unique(dfp$slide), function(si){
  do.call(rbind, lapply(unique(dfp$method), function(mi){
    data.frame(median.value = median(dfp[dfp$slide==si & 
                                           dfp$method==mi,]$value),
               slide = si, method = mi)  
  }))
}))

#-------------------
# get size variables
#-------------------
# make a table with one row per cell
# retain only cell type and size columns, and sample/donor info
dfs <- do.call(rbind, lapply(seq(length(lc)), function(ii){
  message(ii)
  csvi <- lc[[ii]]; csv.fname <- names(lc)[ii]
  cnv <- colnames(csvi)
  # get cell type label
  cnv.type <- colnames(csvi)[grepl(".*Positive$", colnames(csvi))]
  cnv.type <- cnv.type[grepl(paste0(labels, collapse = "|"), cnv.type)]
  do.call(rbind, lapply(cnv.type, function(typei){
    message(typei)
    # get filter terms
    which.type <- which(csvi[,typei]==1) # celltype filter
    csvif <- csvi[which.type,]
    # get size variable
    new.row <- rep("NA", 8)
    names(new.row) <- c(names(cnv.size), "num.cells", "type", "csv.fname")
    for(cni in cnv.size){
      # message(cni)
      cni.filt <- grepl(as.character(cni), colnames(csvif))
      if(length(which(cni.filt))==1){
        datv <- as.numeric(csvif[,cni.filt])
        new.row[names(cnv.size[cnv.size==cni])] <- round(mean(datv),3)
      }
    }
    new.row["num.cells"] <- length(which.type)
    new.row["type"] <- names(labels[labels==gsub("\\..*", "", typei)])
    new.row["csv.fname"] <- csv.fname; new.row
  }))
}))
dfs <- as.data.frame(dfs)
for(c in seq(5)){dfs[,c] <- as.numeric(dfs[,c])}
dfs$slide <- unlist(lapply(dfs$csv.fname, function(fni){
  unlist(strsplit(fni, "_"))[3]
}))
dfs$position <- gsub("[0-9]", "", dfs$slide)
dfs$donor <- gsub("[A-Z]", "", dfs$slide)

# save 
save.fname <- "dfcellsize-byslide_halo.rda"
save(dfs, file = file.path(save.dpath, save.fname))

#-------------
# violin plots
#-------------
dfs <- get(load(file.path(save.dpath, save.fname)))
dfp <- dfs
#colnames(dfp) <- c("akt3", "cell_area", "cyto_area", "nuc_area", 
#                   "nuc_peri", "type", "fname", "donor")
dfp$type <- as.factor(dfp$type)

# get plot data
lgg <- lapply(colnames(dfp)[1:5], function(cni){
  ylab.str <- ifelse(cni=="akt3_expr", "Copies", "Microns")
  dfpi <- dfp[,c(cni, "type")]; dfpi$y = dfpi[,cni]
  ggvp <- ggplot(dfpi, aes(x = type, y = y, color = type)) + 
    geom_violin(draw_quantiles = 0.5) +
    geom_jitter(alpha = 0.5, size = 2) +
    theme_bw() + ylab(ylab.str) +
    ggtitle(cni) + theme(axis.text.x = element_blank(),
                         axis.title.x = element_blank())
  ggvp
})
names(lgg) <- colnames(dfp)[1:5]
# extract legend, and remove from other plots
ggleg <- get_legend(lgg[[1]])
lgg <- lapply(lgg, function(ii){ii+theme(legend.position="none")})

# make multiplot
plot.fname <- "ggvp-agg-bytype_halo-cell-sizes.jpg"
jpeg(file.path(save.dpath, plot.fname), width = 8, height = 4, units = "in",
     res = 400)
grid.arrange(lgg$cell_area, lgg$cytoplasm_area, lgg$nucleus_area,
             lgg$nucleus_perimenter, lgg$akt3_expr, ggleg,
             nrow = 2, bottom = "Cell type", left = "Size")
dev.off()
