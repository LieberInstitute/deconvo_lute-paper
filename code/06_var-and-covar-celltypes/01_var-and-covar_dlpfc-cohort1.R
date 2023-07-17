#!/usr/bin/env R

#
# Variance and covariance matrices, heatmaps before/after cell size factor transform.
#

source("deconvo_method-paper/code/06_var-and-covar-celltypes/00_parameters.R")
sapply(libv, library, character.only = T)

# load mae
mae <- get(load(mae.path))
sce <- mae[[2]]
sce <- logNormCounts(sce)

#---------------------
# get k2 cov, cor data
#---------------------
# get z list, adj and unadj on cell size factor transform
zlist <- z_list(sce, c("glial" = 3, "neuron" = 10), "logcounts", "k2")
# z list by donor
zlist.bydonor <- z_list_bydonor(sce, "Sample", c("glial" = 3, "neuron" = 10), "logcounts", "k2")

# get covariance matrices
zlist.cov <- lapply(zlist, function(item){cov(item)})
zlist.bydonor.cov <- lapply(zlist.bydonor, function(item){
  list("unadj" = cov(item[[1]]), "adj" = cov(item[[2]]))})
names(zlist.bydonor.cov) <- names(zlist.bydonor)

# get correlation matrices
lstat <- covcor_from_zlist(zlist)
lstat.bydonor <- covcor_from_zlist(zlist.bydonor)


Heatmap(lstat$lcov$unadj) + Heatmap(lstat$lcov$adj)

Heatmap(lstat$lcov.diff.unadj.minus.adj)

pheatmap(lstat$lcov.diff.unadj.minus.adj, display_numbers = T)

pheatmap(lstat$lcor.diff.unadj.minus.adj, display_numbers = T)

# distribution of cov difference across donors
cov.diff <- lapply(lstat.bydonor$lcov.diff.unadj.minus.adj, function(item){item[1,2]}) %>% unlist() %>% as.data.frame()
cor.diff <- lapply(lstat.bydonor$lcor.diff.unadj.minus.adj, function(item){item[1,2]}) %>% unlist() %>% as.data.frame()
df.stat.bydonor <- cbind(cov.diff, cor.diff)
colnames(df.stat.bydonor) <- c("cov.diff", "cor.diff")
df.stat.bydonor$sample.id <- rownames(df.stat.bydonor)

# without sample.id labels
ggplot(df.stat.bydonor, aes(x = cov.diff, y = cor.diff)) + 
  theme_bw() + geom_point() + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)

# with sample.id labels
ggplot(df.stat.bydonor, aes(x = cov.diff, y = cor.diff)) + 
  theme_bw() + geom_point() + geom_text(aes(label = sample.id)) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
