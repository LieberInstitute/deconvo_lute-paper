#!/usr/bin/env R

# Author: Sean Maden
# 
# Comparing cell size estimates and size ratios (Non-neuron/Neuron) across HALO
# metrics.
#
#

libv <- c("ggplot2", "ggcorrplot")
sapply(libv, library, character.only = T)

#----------
# load data
#----------
# save filename stem
save.fnstem <- "halo-metrics"

# read output df objects
read.fname <- "dfcellsize_halo.rda"
read.dpath <- save.dpath <- file.path("deconvo_method-paper", "outputs",
                                      "07_cell-size-estimates")
dfc <- get(load(file.path(read.dpath, read.fname)))

#------------
# format data
#------------
# convert to numeric
for(c in seq(5)){dfc[,c] <- as.numeric(dfc[,c])}

# rename cell size variables
colnames(dfc)[seq(5)] <- c("akt3.copies", "cell.area", "cyto.area", "nuc.area", 
                           "nuc.perim")

#------------------------------------------
# correlations by cell type, across samples
#------------------------------------------
ctv <- c("Excit", "Inhib", "Oligo")
lcor <- lapply(ctv, function(cti){
  dfci <- dfc[dfc$type == cti,]
  cor(dfci[,seq(5)], use = "pairwise.complete.obs", method = "spearman")
})
names(lcor) <- ctv

# make correlation heatmaps
for(cti in names(lcor)){
  pdf.fname <- paste0("ggcorrhm-cellsize_",cti,"_",save.fnstem,".pdf")
  pdf.fpath <- file.path(read.dpath, pdf.fname)
  # get plot matrix
  mcori <- lcor[[cti]]
  colnames(mcori) <- gsub(paste0(".", cti), "", colnames(mcori))
  rownames(mcori) <- gsub(paste0(".", cti), "", rownames(mcori))
  # get plot object
  ggcor <- ggcorrplot(mcori, type = "lower", lab = T, title = cti)
  # save new pdf
  ggsave(pdf.fpath, ggcor, width = 5, height = 4,
         device = "pdf", units = "in", dpi = 400)
}

#------------------------------
# compare ratios across samples
#------------------------------
sampv <- unique(dfc$donor)
varv <- colnames(dfc)[2:5]
cnv <- colnames(dfci)[2:5]
dfr <- do.call(rbind, lapply(sampv, function(si){
  dfci <- dfc[dfc$donor==si,]
  # get neuron and glial subsets
  filt.neur <- dfci$type %in% c("Excit", "Inhib")
  dfci.neur <- dfci[filt.neur,]; dfci.gli <- dfci[!filt.neur,]
  # get mean neuron and glial data
  dfri <- unlist(lapply(cnv, function(ci){
    c(mean(dfci.gli[,ci]), mean(dfci.neur[,ci]))
  }))
  names(dfri) <- paste0(rep(cnv, each = 2), ".", c("glial", "neuron"))
  # get ratios
  dfri2 <- data.frame("cell.area.ratio" = dfri[1]/dfri[2],
                      "cyto.area.ratio" = dfri[3]/dfri[4],
                      "nuc.area.ratio" = dfri[5]/dfri[6],
                      "nuc.perim.ratio" = dfri[7]/dfri[8],
                      "donor" = si)
  # return results
  dfr <- cbind(as.data.frame(matrix(dfri, nrow = 1)), dfri2)
  colnames(dfr)[1:8] <- names(dfri)
  return(dfr)
}))

# plot value magnitudes
dfp <- do.call(rbind, lapply(cnv, function(ci){
  filt.cnv <- grepl(paste0(ci, ".ratio"), colnames(dfr))
  dfi <- dfr[,filt.cnv]
}))
ggplot(dfr, aes(x = ))

# get correlation matrix
mcor <- cor(dfr[,seq(4)], use = "pairwise.complete.obs", method = "spearman")

# make new corr figure
pdf.fname <- paste0("ggcorrhm-cellratio_",save.fnstem,".pdf")
pdf.fpath <- file.path(read.dpath, pdf.fname)
# get plot matrix
mcori <- lcor[[cti]]
colnames(mcori) <- gsub(paste0(".", cti), "", colnames(mcori))
rownames(mcori) <- gsub(paste0(".", cti), "", rownames(mcori))
# get plot object
ggcor <- ggcorrplot(mcori, type = "lower", lab = T, 
                    title = "Ratio (glial/neuron)")
# save new pdf
ggsave(pdf.fpath, ggcor, width = 5, height = 4,
       device = "pdf", units = "in", dpi = 400)
