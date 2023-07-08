#!/usr/bin/env R

#
#
#
#

# load data
source("deconvo_method-paper/code/10_donor-bias-simulations-continued/00_parameters.R")
sapply(libv, library, character.only = T)
sce <- get(load(sce.path))

# means and variances by donor
sce[[new.group.variable]] <- group.vector <- 
  paste0(sce[[k.marker.variable]], ";", sce[[sample.variable.name]])

unique.groups <- group.vector %>% unique() # unique(sce[[new.group.variable]])

# variance by types

expression.matrix <- sce[,sce.filter] %>% assays()[[assay.name]]

variance.matrix <- aggregate()

sce.var <- do.call(cbind, lapply(unique.groups, function(groupi){
  message("working on group level ", groupi)
  sce.filter <- sce[[varname]] == groupi
  rowVars(assays(sce[,sce.filter])[[assayname]])
}))

sce.mean <- do.call(cbind, lapply(unique.groups, function(groupi){
  message("working on group level ", groupi)
  filt <- scei[[varname]] == groupi
  rowMeans(assays(scei[,filt])[[assayname]])
}))

colnames(sce.mean) <- colnames(sce.var) <- unique.groups

# save
# variances
fname <- paste0("matrix-donor-type_variances_",markeri, "_", handle.str,".rda")
save(sce.var, file = file.path(save.dpath, fname))
# means
fname <- paste0("matrix-donor-type_means_",markeri, "_", handle.str,".rda")
save(sce.mean, file = file.path(save.dpath, fname))

#-------------------
# plot mean-variance
#-------------------
genev <- rownames(sce.mean)
groupv <- colnames(sce.mean)
dfp <- do.call(rbind, lapply(seq(ncol(sce.mean)), function(ii){
  dfpi <- data.frame(gene = genev, mean = sce.mean[,ii], var = sce.var[,ii])
  dfpi$group <- groupv[ii]; return(dfpi)
}))
dfp$celltype <- gsub(";.*", "", dfp$group)
dfp$donor <- gsub(".*;", "", dfp$group)

ggsm <- ggplot(dfp, aes(x = mean, y = var, color = donor, group = donor)) +
  geom_smooth() + facet_wrap(~celltype) + theme_bw() +
  geom_abline(intercept = 0, slope = 1, color = "black") +
  scale_x_log10() + scale_y_log10()

# save jpg
save.dpath <- "."
fname <- "ggsmooth-composite_k2-types-by-donor_ro1-dlpfc.jpg"
jpeg(file = file.path(save.dpath, fname), width = 8, height = 4, 
     units = "in", res = 400)
print(ggsm); dev.off()

#---------------------------------
# get random background dispersion
#---------------------------------


