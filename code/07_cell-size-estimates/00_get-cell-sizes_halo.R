#!/usr/bin/env R

# Get the cell sizes from HALO outputs.
#
#

library(ggplot2)

#----------
# load data
#----------
# manage paths
base.path <- "."
save.dname <- "07_cell-size-estimates"
save.dpath <- file.path(base.path, "deconvo_method-paper", "outputs", save.dname)

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

#-----------
# set params
#-----------
labels = c("Endo" = "CLDN5", "Astro" = "GFAP",
           "Inhib" = "GAD1", "Excit" = "SLC17A7",
           "Micro" = "TMEM119", "Oligo" = "OLIG2")

cnv.size <- c("Cell.Area.*", "Cytoplasm.Area.*", "Nucleus.Area.*",
              "Nucleus.Perimeter.*", "AKT3..Opal.570..Copies")

# get flattened, filtered list
lc1 <- lapply(lcsv[[1]], function(ii){ii})
lc2 <- lapply(lcsv[[2]], function(ii){ii})
lc <- c(lc1[grepl("Final", names(lc1))], lc2[grepl("Final", names(lc2))])

#-------------------
# get size variables
#-------------------
dfsize <- do.call(rbind, lapply(seq(length(lc)), function(ii){
  message(ii)
  csvi <- lc[[ii]]; csv.fname <- names(lc)[ii]
  cnv <- colnames(csvi)
  # get cell type label
  cnv.type <- colnames(csvi)[grepl(".*Positive$", colnames(csvi))]
  cnv.type <- cnv.type[grepl(paste0(labels, collapse = "|"), cnv.type)]
  do.call(rbind, lapply(cnv.type, function(typei){
    message(typei)
    which.type <- csvi[,typei]==1
    cnvf <- cnv[grepl(paste0(cnv.size,collapse ="|"), cnv)]
    cnv.size.in <- cnv.size[cnv.size %in% cnvf]
    cnv.size.out <- cnv.size[!cnv.size %in% cnvf]
    # parse available colnames
    datv <- csvi[which.type, cnv.size.in, drop = F]
    if(ncol(datv) > 1){datv <- colMeans(datv)}
    # parse unavailable colnames
    for(c in cnv.size.out){datv[c] <- "NA"}
    # use standard order
    datv <- datv[order(names(datv), cnv.size)]
    datv["type"] <- names(labels[labels==gsub("\\..*", "", typei)])
    datv["csv.fname"] <- csv.fname; datv
  }))
}))
for(c in seq(5)){dfsize[,c] <- as.numeric(dfsize[,c])}
dfsize$donor <- unlist(lapply(dfsize$csv.fname, function(fni){
  unlist(strsplit(fni, "_"))[3]
}))

# save 
save.fname <- "dfcellsize_halo.rda"
save(dfsize, file = file.path(save.dpath, save.fname))

#-------------
# violin plots
#-------------
dfp <- dfsize
colnames(dfp) <- c("akt3", "cell_area", "cyto_area", "nuc_area", 
                   "nuc_peri", "type", "fname", "donor")
dfp$type <- as.factor(dfp$type)

lgg <- lapply(colnames(dfp)[1:5], function(cni){
  dfpi <- dfp[,c(cni, "type")]; dfpi$y = dfpi[,cni]
  ggvp <- ggplot(dfpi, aes(x = type, y = y, color = type)) + 
    geom_violin(draw_quantiles = 0.5) +
    theme_bw() + ylab(cni)
  ggvp
})
names(lgg) <- colnames(dfp)[1:5]
