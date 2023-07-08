
library(lute)

outputs.path <- here("deconvo_method-paper", "outputs", "18_updated-marker-sets-analyses")

# using concordant, overlapping markers
sce.name <- "sce_top-concordant-overlapping_training.rda"
sce.path <- here(outputs.path, sce.name)
sce <- get(load(sce.path))
results <- deconvolution.experiment(sce = sce,
                                    s = c("glial" = 3, "neuron" = 10),
                                    assay.name = "counts",
                                    sample.id.variable = "Sample",
                                    celltype.variable = "k2")

# using old approach markers
se.name <- "se-adj_marker-filter-all_train.rda"
se.path <- here(outputs.path, se.name)
se <- se.train.old <- get(load(se.path))
sce <- SingleCellExperiment(se)
colData(sce) <- colData(se)
assays(sce) <- list(counts = as.matrix(assays(se)[["counts_adj"]]))
sce <- scuttle::logNormCounts(sce)
results <- deconvolution.experiment(sce = sce,
                                    s = c("glial" = 3, "neuron" = 10),
                                    assay.name = "counts",
                                    sample.id.variable = "Sample",
                                    celltype.variable = "k2")

results$plots.list$neuron$proportions.scatterplot
results$plots.list$neuron$abs.error.jitterbox

#---------------
# results in mrb
#---------------
sce.mrb.name <- "sce-mrb_dlpfc.rda"
sce.mrb.path <- here("deconvo_method-paper", "outputs", "09_manuscript", sce.mrb.name)
sce.mrb <- get(load(sce.mrb.path))

markers.old.top40 <- metadata(se.train.old)$top.40.markers.results$markers
sce.mrb.old <- sce.mrb[rownames(sce.mrb) %in% markers.old.top40,]
celltype.variable.vector <- ifelse(grepl("Excit|Inhib", sce[["cellType"]]),"neuron",
                                   ifelse(sce[["cellType"]] %in% 
                                            c("Oligo", "OPC", "Micro", "Astro"), 
                                          "glial", "other"))
sce.mrb.old[["k2"]] <- celltype.variable.vector
filter.sce <- !sce.mrb.old[["k2"]] == "other"
sce.mrb.old.filter <- sce.mrb.old[,filter.sce]

results.mrb.old <- deconvolution.experiment(sce = sce.mrb.old.filter,
                                    s = c("glial" = 3, "neuron" = 10),
                                    assay.name = "counts",
                                    sample.id.variable = "donor",
                                    celltype.variable = "k2")


results.mrb.old <- deconvolution.experiment(sce = sce.mrb.old.filter,
                                            s = c("glial" = 3, "neuron" = 10),
                                            assay.name = "counts",
                                            sample.id.variable = "donor",
                                            celltype.variable = "k2")



