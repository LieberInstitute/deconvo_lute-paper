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
read.fname <- ""
read.dpath <- out.dpath <- file.path("deconvo_method-paper", "outputs", 
                                     "07_cell-size-estimates")
dfm <- get(load(file.path(read.dpath, read.fname)))

#------------
# format data
#------------
# convert to numeric
for(c in seq(36)){dfm[,c] <- as.numeric(dfm[,c])}

#------------------------------------------
# correlations by cell type, across samples
#------------------------------------------
ctv <- c("Excit", "Inhib", "Oligo")
lcor <- lapply(ctv, function(cti){
  dffm <- dfm[,grepl(cti, colnames(dfm))]
  for(c in seq(ncol(dffm))){dffm[,c] <- as.numeric(dffm[,c])}
  cor(dffm, use = "pairwise.complete.obs", method = "spearman")
})
names(lcor) <- ctv

# make correlation heatmaps
for(cti in names(lcor)){
  pdf.fname <- paste0("ggcorrhm-cellsize_",cti,
                      "_",save.fnstem,".pdf")
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

#-----------------------------
# compare means across samples
#-----------------------------
ctv <- c("Excit", "Inhib", "Oligo")
meanv <- colMeans(dfm[,seq(36)], na.rm = T)
dfp <- do.call(rbind, lapply(ctv, function(cti){
  mfv <- meanv[grepl(cti, names(meanv))]
  mfv.halo <- mfv[grepl("^halo\\..*", names(mfv))]
  dfp <- as.data.frame(matrix(mfv.halo, ncol = 1))
  dfp$halo.var <- gsub(paste0("\\.", cti, "$"), "", names(mfv.halo))
  dfp$celltype <- cti
  dfp$sce.mean.across.samples <- mfv[grepl("^sce\\..*", names(mfv))]
  return(dfp)
}))
colnames(dfp)[1] <- "halo.mean.across.samples"

# get ggplot
# make scatterplot
ggpt <- ggplot(dfp, aes(x = sce.mean.across.samples, 
                        y = halo.mean.across.samples,
                        color = celltype, shape = halo.var)) + 
  geom_point(size = 5, alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, col = "black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("SCE") + ylab("HALO") + ggtitle("Mean across samples")
# save scatterplot
pdf.fname <- paste0("ggpt-cellsize_", save.fnstem, ".pdf")
pdf.fpath <- file.path(read.dpath, pdf.fname)
ggsave(pdf.fpath, ggpt, width = 5, height = 4, 
       device = "pdf", units = "in", dpi = 400)

# make facet plot
ggpt <- ggplot(dfp, aes(x = sce.mean.across.samples, 
                        y = halo.mean.across.samples,
                        color = celltype)) + geom_point() +
  geom_abline(intercept = 0, slope = 1, col = "black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("SCE") + ylab("HALO") + ggtitle("Mean across samples")
ggf <- ggpt + facet_wrap(~halo.var)
# save facet plot
pdf.fname <- paste0("ggptfacet-cellsize_", save.fnstem, ".pdf")
pdf.fpath <- file.path(read.dpath, pdf.fname)
ggsave(pdf.fpath, ggf, width = 6.5, height = 3.5, 
       device = "pdf", units = "in", dpi = 400)


