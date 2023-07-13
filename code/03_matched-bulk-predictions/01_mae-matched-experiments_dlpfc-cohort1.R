#!/usr/bin/env R

#
#
#

source("deconvo_method-paper/code/03_matched-bulk-predictions/00_parameters-matched-assays.R")
sapply(libv, library, character.only = T)

# load mae 
mae.filepath <- here("deconvo_method-paper", "outputs", "01_prepare-datasets", "mae_final.rda")
mae <- get(load(mae.filepath))

#-----------------------------
# get cell size factor series
#-----------------------------
# load cell size scale factors
df.csf <- get_csf_reference()
# get cell size scale factors from mae data

#--------------
# k2 experiment
#--------------
# get_ymatch_experiment_results()

sample.id.vector <- unique(colData(mae)[,1])
sample.id <- sample.id.vector[1]

s.vector <- c("glial" = 2, "neuron" = 3)
sce.iter <- mae[[2]]

rownames(rse.all) <- rowData(rse.all)$Symbol
filter.y.sample <- colData(rse.all)[,"batch.id2"]==sample.id
filter.y.marker <- rownames(rse.all) %in% rownames(sce.iter)
rse.iter <- rse.all[filter.y.marker, filter.y.sample]

s.vector.pred <- order_svector(s.vector.pred)

# get results for a single iteration of an experiment
# use with get_ypb_experiment_series()
if(assay.name == "logcounts" & !"logcounts" %in% names(assays(sce))){sce <- scuttle::logNormCounts(sce)}
unique.sample.id.vector <- unique(sce[[sample.id.variable]])
dfp <- do.call(rbind, lapply(unique.sample.id.vector, function(sample.id){
  sce.iter <- sce[,sce[[sample.id.variable]]==sample.id]
  y.iter <- mae[[1]] %>% as.matrix()
  prop.true.iter <- table(sce.iter[[celltype.variable]]) %>% prop.table() %>% as.matrix() %>% t()
  prop.pred.iter <- lute(sce = sce.iter, y = ypb.iter, assay.name = assay.name, 
                         celltype.variable = celltype.variable, s = s.vector.pred, 
                         typemarker.algorithm = NULL, return.info = FALSE,
                         deconvolution.algorithm = deconvolution.algorithm)$deconvolution.results
  colnames(prop.pred.iter) <- paste0(colnames(prop.pred.iter), ".pred")
  colnames(prop.true.iter) <- paste0(colnames(prop.true.iter), ".true")
  dfp.iter <- cbind(prop.true.iter, prop.pred.iter) %>% as.data.frame()
  Sys.sleep(system.sleep.sec)
  dfp.iter
}))
rownames(dfp) <- unique.sample.id.vector





assay.name <- "logcounts"

s.vector.pred()

lute(sce = sce.iter, y = ypb.iter, assay.name = assay.name, 
     celltype.variable = celltype.variable, s = s.vector.pred, 
     typemarker.algorithm = NULL, return.info = FALSE,
     deconvolution.algorithm = deconvolution.algorithm)$deconvolution.results


# get experiment results tables
dfp.tall <- get_ypb_experiment_series(sce.k2, sample.id.variable = "Sample", 
                                      celltype.variable = "k2", assay.name = "logcounts",
                                      s.vector = c("glial" = 3, "neuron" = 10),
                                      algorithm.name = "nnls", return.dimensions = "tall")
