#!/usr/bin/env R

# Author: Sean Maden
#
# Test independent pseudobulk from DLPFC cohort1 snRNAseq data.
#
#
#
#
#

source("scripts/02_pseudobulk/00_param.R")
sapply(libv, library, character.only = T)
list.sce.markers <- get(load(sce.markers.list.path))
sce <- list.sce.markers$k4
s.vector <- c("Excit" = 10, "Inhib" = 10, "non_oligo_glial" = 3, "Oligo" = 3)

# get experiment results tables
unique.sample.id.vector <- unique(sce[["Sample"]])
sample.id.filter1 <- sce[["Sample"]] %in% unique.sample.id.vector[1]
sample.id.filter2 <- sce[["Sample"]] %in% unique.sample.id.vector[2]
sample.id.filter3 <- sce[["Sample"]] %in% unique.sample.id.vector[3]
sample.id.filter4 <- sce[["Sample"]] %in% unique.sample.id.vector[4]
sample.id.filter5 <- sce[["Sample"]] %in% unique.sample.id.vector[5]

# get dfp.tall with memory management
dfp.tall1 <- get_ypb_experiment_series(sce[,sample.id.filter1], sample.id.variable = "Sample", 
                                       celltype.variable = "k4", assay.name = "logcounts",
                                       s.vector = s.vector,
                                       algorithm.name = "nnls", return.dimensions = "tall",
                                       system.sleep.sec = 10)
dfp.tall2 <- get_ypb_experiment_series(sce[,sample.id.filter2], sample.id.variable = "Sample", 
                                       celltype.variable = "k4", assay.name = "logcounts",
                                       s.vector = s.vector,
                                       algorithm.name = "nnls", return.dimensions = "tall",
                                       system.sleep.sec = 10)
dfp.tall3 <- get_ypb_experiment_series(sce[,sample.id.filter3], sample.id.variable = "Sample", 
                                       celltype.variable = "k4", assay.name = "logcounts",
                                       s.vector = s.vector,
                                       algorithm.name = "nnls", return.dimensions = "tall",
                                       system.sleep.sec = 10)
dfp.tall4 <- get_ypb_experiment_series(sce[,sample.id.filter4], sample.id.variable = "Sample", 
                                       celltype.variable = "k4", assay.name = "logcounts",
                                       s.vector = s.vector,
                                       algorithm.name = "nnls", return.dimensions = "tall",
                                       system.sleep.sec = 10)
dfp.tall5 <- get_ypb_experiment_series(sce[,sample.id.filter5], sample.id.variable = "Sample", 
                                       celltype.variable = "k4", assay.name = "logcounts",
                                       s.vector = s.vector,
                                       algorithm.name = "nnls", return.dimensions = "tall",
                                       system.sleep.sec = 10)


list.dfp.tall <- lapply(unique.sample.id.vector, function(sample.id){
  
  filter.sample.id <- sce$Sample==sample.id
  
  dfp.tall.iter <- get_ypb_experiment_series(sce[,filter.sample.id], 
                                             sample.id.variable = "Sample",
                                             celltype.variable = "k4", 
                                             assay.name = "logcounts",
                                             s.vector = s.vector,
                                             algorithm.name = "nnls", 
                                             return.dimensions = "tall",
                                             system.sleep.sec = 5)
  
  return(dfp.tall.iter)
  
})

unique.sample.id.vector <- unique(sce[["Sample"]])

sample.id <- unique.sample.id.vector[1:2]

filter.sample.id <- sce$Sample %in% sample.id

dfp.tall1 <- get_ypb_experiment_series(sce[,filter.sample.id], sample.id.variable = "Sample", 
                                       celltype.variable = "k4", assay.name = "logcounts",
                                       s.vector = s.vector,
                                       algorithm.name = "nnls", return.dimensions = "tall",
                                       system.sleep.sec = 5)

gc()


dfp.tall1 <- get_ypb_experiment_series(sce[,sample.id.filter1], sample.id.variable = "Sample", 
                                       celltype.variable = "k4", assay.name = "logcounts",
                                       s.vector = s.vector,
                                       algorithm.name = "nnls", return.dimensions = "tall",
                                       system.sleep.sec = 5)
gc()

dfp.tall2 <- get_ypb_experiment_series(sce[,sample.id.filter2], sample.id.variable = "Sample", 
                                       celltype.variable = "k4", assay.name = "logcounts",
                                       s.vector = s.vector,
                                       algorithm.name = "nnls", return.dimensions = "tall",
                                       system.sleep.sec = 5)
gc()

dfp.tall3 <- get_ypb_experiment_series(sce[,sample.id.filter3], sample.id.variable = "Sample", 
                                       celltype.variable = "k4", assay.name = "logcounts",
                                       s.vector = s.vector,
                                       algorithm.name = "nnls", return.dimensions = "tall",
                                       system.sleep.sec = 5)
gc()
dfp.tall <- rbind(dfp.tall1, rbind(dfp.tall2, dfp.tall3))

# get dfp.wide with memory management
dfp.wide1 <- get_ypb_experiment_series(sce[,sample.id.filter1], sample.id.variable = "Sample", 
                                      celltype.variable = "k4", assay.name = "logcounts",
                                      s.vector = s.vector,
                                      algorithm.name = "nnls", return.dimensions = "wide")
gc()
dfp.wide2 <- get_ypb_experiment_series(sce[,sample.id.filter2], sample.id.variable = "Sample", 
                                       celltype.variable = "k4", assay.name = "logcounts",
                                       s.vector = s.vector,
                                       algorithm.name = "nnls", return.dimensions = "wide")
gc()
dfp.wide3 <- get_ypb_experiment_series(sce[,sample.id.filter2], sample.id.variable = "Sample", 
                                       celltype.variable = "k4", assay.name = "logcounts",
                                       s.vector = s.vector,
                                       algorithm.name = "nnls", return.dimensions = "wide")
gc()
dfp.wide <- rbind(dfp.wide1, rbind(dfp.wide2, dfp.wide3))

dfp.ct <- dfp_tall_by_celltype(dfp.wide) # dfp.wide, tall by cell type

# make new plots
# plot proportions multipanel
ggplot(dfp.ct, aes(x = true.noscale, y = pred.noscale)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_hline(yintercept = 0.5) + geom_vline(xintercept = 0.5) + theme_bw() +
  xlab("True proportion") + ylab("Predicted proportion") +
  xlim(0, 1) + ylim(0, 1) + facet_wrap(~celltype)

# jitterbox -- jittered points and boxplots of absolute errors
dfp.ae1 <- data.frame(celltype = dfp.ct$celltype,
                      abs.error = dfp.ct$abs.error.withscale)
dfp.ae1$type <- "withscale"
dfp.ae2 <- data.frame(celltype = dfp.ct$celltype,
                      abs.error = dfp.ct$abs.error.noscale)
dfp.ae2$type <- "noscale"
dfp.ae <- rbind(dfp.ae1, dfp.ae2)
ggplot(dfp.ae, aes(x = celltype, y = abs.error)) + geom_jitter(alpha = 0.5) + 
  geom_boxplot(color = "cyan", alpha = 0) + theme_bw() + facet_wrap(~type)

# save environment
save.image("./env/02_pseudobulk/03_k4.RData")

