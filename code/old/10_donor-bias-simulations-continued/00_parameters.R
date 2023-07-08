#-------
# header
#-------
libv <- c("lute", "here", "dplyr", "SingleCellExperiment", "limma", 
    "SummarizedExperiment", "scuttle", "ggplot2", "splatter", "scater", 
    "gridExtra", "ggpubr", "GGally")
sapply(libv, library, character.only = TRUE)

# helper functions
get_sce_donor_bias <- function(donor.bias.coeff = 1, mean.vector = NULL,
                               total.cells = 10, num.genes.iter = 10,
                               num.types = 1, seed.num = 0){
  set.seed(seed.num)
  if(is(mean.vector, "NULL")){mean.vector = seq(1, 30, 1)}
  sce <- do.call(rbind, lapply(mean.vector, function(meani){
    d <- donor.bias.coeff*(meani/(meani^2))
    scei <- random_sce(num.types = num.types, num.cells = total.cells, 
                       expr.mean = meani, num.genes = num.genes.iter, 
                       dispersion = d)
    return(scei)
  }))
  return(as(sce, "SingleCellExperiment"))
}

simulate_sce_donor_bias <- function(donor.coeff.vector = NULL,
                                    mean.vector = NULL, 
                                    cells.per.donor = 1000,
                                    num.genes.iter = 10){
  lr <- list()
  if(is(donor.coeff.vector, "NULL")){
    donor.coeff.vector <- rnorm(10, mean = 200, sd = 30)
  }
  sce <- do.call(cbind, lapply(donor.coeff.vector, function(ii){
    scei <- get_sce_donor_bias(donor.bias.coeff = ii, 
                               mean.vector = mean.vector, 
                               total.cells = cells.per.donor,
                               num.genes.iter = num.genes.iter)
    scei[["donor.coeff"]] <- ii; scei
  }))
  dfp <- do.call(rbind, lapply(donor.coeff.vector, function(ii){
    ctf <- counts(sce[,sce[["donor.coeff"]]==ii])
    dfi <- data.frame(mean = rowMeans(ctf), var = rowVars(ctf))
    dfi$donor.coeff <- ii; dfi
  }))
  dfp$donor.coeff <- as.character(dfp$donor.coeff)
  ggsm <- ggplot(dfp, aes(x = mean, y = var, color = donor.coeff)) + 
    geom_abline(slope = 1, intercept = 0) + 
    scale_x_log10() + scale_y_log10() + geom_smooth() +
    theme(legend.position = 'none') +
    xlab('Mean (log10)') + ylab("Var (log10)")
  lr[["sce"]] <- sce; lr[["dfp"]] <- dfp; lr[["ggsmooth"]] <- ggsm
  return(lr)
}

donors_from_sce <- function(sce, bias.vector = seq(3), bias.means = seq(3)){
  ct <- counts(sce)
  medians <- rowMedians(ct) + 1
  num.cells <- ncol(ct)
  num.genes <- nrow(ct)
  sce.new <- do.call(cbind, lapply(bias.vector, function(biasi){
    sce.new <- sce
    d <- biasi*(medians/(medians^2))
    ct.new <- do.call(rbind, lapply(seq(num.genes), function(ii){
      rnbinom(num.cells, size = d[ii], mu = medians[ii])
    }))
    rownames(ct.new) <- rownames(ct)
    colnames(ct.new) <- colnames(ct)
    assays(sce.new)[["counts"]] <- ct.new
    sce.new[["donor.bias"]] <- biasi
    return(sce.new)
  }))
  return(sce.new)
}

plot_mean_var_sce <- function(sce, donor.variable = NULL){
  if(is(donor.variable, "NULL")){
    ct <- counts(sce)
    dfp <- data.frame(mean = rowMeans(ct), var = rowVars(ct))
    ggplot(dfp, aes(x = mean, y = var)) +
      geom_smooth() + geom_abline(intercept = 0, slope = 1) +
      theme_bw() + scale_x_log10() + scale_y_log10()
  } else{
    donor.vector <- unique(sce[[donor.variable]])
    dfp <- do.call(rbind, lapply(donor.vector, function(donori){
      ctf <- counts(sce[,sce[[donor.variable]]==donori])
      dfi <- data.frame(mean = rowMeans(ctf), var = rowVars(ctf))
      dfi$donor <- donori; dfi
    }))
    dfp$donor <- as.character(dfp$donor)
    ggplot(dfp, aes(x = mean, y = var, color = donor)) +
      geom_smooth() + geom_abline(intercept = 0, slope = 1) +
      theme_bw() + scale_x_log10() + scale_y_log10()
  }
}

# load
project.handle.string <- "ro1-dlpfc"
sce.list.file.name <- paste0("list-scef_markers-k2-k3-k4_",project.handle.string,".rda")
sce.list.path <- file.path(sce.load.path, sce.list.file.path)

workflow.table.filename <- paste0("workflow-table-all_", save.filename.stem, ".csv")
workflow.table.path <- here(base.path, workflow.table.filename)

#-----------------------
# parameters for scripts
#-----------------------
assay.name <- "count"

# mean-variance (01)
sample.variable.name <- "Sample"
k.marker.variable <- "k2"
experiment.name <- "ro1-dlpfc"
assayname <- "counts_adj"
new.group.variable <- varname <- "group.variable.name"
# sce load path
sce.load.path <- here("deconvo_method-paper", "outputs", "09_manuscript")
sce.file.name <- paste0("sce_marker-adjustment-",k.marker.variable,"_",experiment.name,".rda")
sce.path <- here(sce.load.path, sce.file.name)

# donor bias dispersion (02)
donor.variable.name <- "donor.bias"
seed.number <- 1

# cells by sample (03)
celltype.variable <- "cellType"

# subsample visualizations (04)
# get save path
code.dname <- "10_donor-bias-simulations-continued"
proj.dname <- "deconvo_method-paper"
save.dpath <- file.path(proj.dname, "outputs", code.dname)
# get base path
rnf.dname <- "r-nf_deconvolution_intra-sample-bias"
base.path <- file.path(save.dpath, rnf.dname, "data")
# read tables
stemv <- paste0(c("inter", "intra"), "-sample-bias")
base.stem <- "r-nf_deconvolution_"
dnv <- paste0(base.stem, stemv)
# experiment design groups
across.batch.label <- "between-batch"
within.batch.label <- "within-batch"

# across-donor bias, donor-bias simulations (04)
# experiment parameters
# get load path
load.code.directory.name <- "09_manuscript"
project.directory.name <- "deconvo_method-paper"
load.dpath <- file.path(project.directory.name, "outputs", code.directory.name)
# get save path
save.code.directory.name <- "10_donor-bias-simulations-continued"
save.directory.path <- file.path(project.directory.name, "outputs", code.directory.name)
# variable values for experiment
seed.num <- 0
iterations <- 1000
fraction.cells <- 25
num.sample.iter <- 3
count.min <- 200
num.batch <- 200
celltype.variable <- k.marker.variable <- "k2"
group.variable <- "Sample"
assay.name <- "counts_adj"
save.fnstem <- paste0("inter-sample_", project.handle.string)
save.filename.stem <- paste0("inter-sample_", experiment.name)
scale.factor <- c("glial" = 3, "neuron" = 10)
methodv <- c("nnls", "music", "epic", "deconrnaseq")
base.directory.workflow.data <- "data"
base.path.workflow.data <- file.path(save.dpath, rnf.dname, base.directory.workflow.data)
rnf.workflow.directory <- "r-nf_deconvolution_inter-sample-bias"
workflow.path <- here(save.directory.path, rnf.workflow.directory)

# manuscript figures, across-donor bias, donor-bias simulations (05)
results.table.filter.string <- "results_table_.*"
data.dpath <- base.path.workflow.data
results.table.file.path <- file.path(save.dpath, rnf.dname, "data", rt.fname)
variable.column.names <- c("prop.pred.type1", "prop.pred.type2", "bias.type1", "bias.type2", "rmse.types")
iterations.index.column.name <- "iterations_index"
algorithm.column.name <- "deconvolution_method"
jitter.plot.xaxis.label <- "Deconvolution.algorithm"
# plot parameters
violin.plots.five.filename <- "ggvp-comp_inter-sample-bias.jpg"
violin.plots.two.filename <- "ggvp-comp-bias_inter-sample-bias.jpg"
scatterplot.filename <- "ggpt-facet-method_inter-sample-bias.jpg"

# within-batch bias simulations (06)
# get load path
code.dname <- "09_manuscript"
proj.dname <- "deconvo_method-paper"
load.dpath <- file.path(proj.dname, "outputs", code.dname)
# get save path
code.dname <- "10_donor-bias-simulations-continued"
proj.dname <- "deconvo_method-paper"
save.dpath <- file.path(proj.dname, "outputs", code.dname)
# get params for experiment
seed.num <- 0
celltype.variable <- "k2"
proj.handle <- "ro1-dlpfc"
save.fnstem <- paste0("intra-sample_", proj.handle)
group.variable <- "Sample"
assay.name <- "counts_adj"
methodv <- c("nnls", "music", "epic", "deconrnaseq")
iterations <- 1000
# fraction.cells <- 25
num.sample.iter <- 3
scale.factor <- c("glial" = 3, "neuron" = 10)
results.table.withingroup.filename.stem <- "r-nf_deconvolution_intra-sample-bias"
base.path <- "data"
base.path <- file.path(save.dpath, rnf.dname, base.path)
# save data
which.save = c("li", "sce")
save.names = list(sce.name = "sce.rda", li.name = "lindex.rda")
# load data
fname <- paste0("list-scef_markers-k2-k3-k4_",proj.handle,".rda")
fpath <- file.path(load.dpath, fname)
count.minimum.acrossgroup <- 200
rnf.workflow.directory <- "r-nf_deconvolution_inter-sample-bias"
workflow.path.withingroup <- here(save.directory.path, rnf.workflow.directory)
save.filename.stem.withingroup <- paste0("intra-sample_", proj.handle)
wt.fnamei <- paste0("workflow-table-all_",save.filename.stem.withingroup,".csv")
workflow.table.path.withingroup <- file.path(base.path, wt.fnamei)
num.batch.withingroup <- 200

# plot within-batch bias (07)
results.filt <- "results_table_.*"
data.dpath <- file.path(save.dpath, rnf.dname, "data")