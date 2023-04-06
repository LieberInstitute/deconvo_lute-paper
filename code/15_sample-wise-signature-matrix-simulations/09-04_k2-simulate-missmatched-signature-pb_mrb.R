
# setup simulation
# perform deconvolution
# uses missmatched pseudobulk and signature matrices within samples
# includes null cell sizes adjustment

source("deconvo_method-paper/code/15_sample-wise-signature-matrix-simulations/00_parameters-script-set-15.R")
sapply(libv, library, character.only = T)
sce.path <- here("deconvo_method-paper", "outputs", "09_manuscript", "sce-mrb_dlpfc.rda")
sce <- get(load(sce.path))

# get k2 labels and filter
cell.type.vector <- sce[["cellType"]]
k2.cell.type.vector <- ifelse(grepl("Excit|Inhib", cell.type.vector), "neuron",
                              ifelse(cell.type.vector %in% c("Oligo", "OPC", "Astro", "Micro"), 
                                     "glial", "other"))
sce[["k2"]] <- k2.cell.type.vector
filter.other <- !sce[["k2"]]=="other"
sce <- sce[,filter.other]

# filter markers
reference.list.path <- here("deconvo_method-paper", "outputs", 
                            "15_k2-simulations_within-sample-matched", 
                            "list-references-by-overlap-rate_dlpfc-ro1-train.rda")
reference.list <- get(load(reference.list.path))
markers.vector <- rownames(reference.list[[3]][[1]])
sce <- sce[markers.vector,]

# experiment parameters
sample.id.variable <- "donor"
celltype.variable <- "k2"
assay.name <- "counts"
s <- c("glial" = 3, "neuron" = 10)
deconvolution.algorithm <- "nnls"
summary.method  <- "mean"

experiment <- deconvolution.experiment.permute.groups(
  sample.id.variable = sample.id.variable, 
  celltype.variable = celltype.variable, assay.name = assay.name, s = s,
  deconvolution.algorithm = deconvolution.algorithm, 
  summary.method = summary.method, verbose = TRUE)

deconvolution.experiment.permute.groups <- function(verbose = FALSE, ...){
  if(verbose){message("running permutation experiments...")}
  group.id.vector <- sce[[sample.id.variable]]
  unique.group.id.vector <- group.id.vector %>% unique() %>% as.character()
  results.table.list <- lapply(unique.group.id.vector, function(group.id.z){
    if(verbose){message("working on group id ", group.id.z, "...")}
    if(verbose){message("getting main signature matrix...")}
    filter.sce <- sce[[sample.id.variable]]==group.id.z
    sce.group.z <- sce[,filter.sce]
    z <- signature_matrix_from_sce(sce.group.z,
                                   celltype.variable = celltype.variable,
                                   summary.method = summary.method,
                                   assay.name = assay.name)
    
    if(verbose){message("getting pseudobulks...")}
    pb.filter <- !unique.group.id.vector==group.id.z
    unique.group.id.pb <- unique.group.id.vector[pb.filter] %>% unique()
    ypb.list <- lapply(unique.group.id.pb, function(group.id){
      if(verbose){
        message("working on pseudobulk for sample id ", group.id, "...")}
      filter.group <- sce[[sample.id.variable]]==group.id
      ypb_from_sce(sce = sce[,filter.group], assay.name = assay.name, 
                   celltype.variable = celltype.variable, S = s)
    })
    ypb.table <- do.call(cbind, ypb.list) %>% as.data.frame()
    colnames(ypb.table) <- unique.group.id.pb
    
    if(verbose){message("getting experiment series...")}
    experiment <- deconvolution.experiment(sce = sce, 
                                        y = ypb.table, s = s, z = z, 
                                        assay.name = assay.name, 
                                        sample.id.variable = sample.id.variable, 
                                        experiment.labels = group.id.z,
                              deconvolution.algorithm = deconvolution.algorithm,
                              celltype.variable = celltype.variable)
    
    if(verbose){message("finished with group id, returning results table...")}
    results.table.iteration <- experiment$results.table
    results.table.iteration$group.id.signature <- group.id.z
    results.table.iteration
  })
  results.table <- do.call(rbind, results.table.list) %>% as.data.frame()
  # make plots
  plots.list <- deconvolution.results.plots.permutations(results.table)
  return(
    list(
      results.table = results.table, plots.list = plots.list))
}

deconvolution.results.plots.permutations <- function(results.table, verbose = FALSE){
  if(verbose){"Making results plots for permutation experiments..."}
  filter.results <- results.table$sample.id==results.table$group.id.signature
  filter.types.vector <- c("matched", "miss-matched", "all")
  lgg.permute <- lapply(filter.types.vector, function(filter.type){
    if(filter.type=="matched"){
      filter.final <- filter.results
    } else if(filter.type=="miss-matched"){
      filter.final <- !filter.results
    } else{
      filter.final <- seq(nrow(results.table))
    }
    results.filtered <- results.table[filter.final,]
    results.filtered$sample.id <- results.filtered$group.id.signature
    lgg.filter <- deconvolution.results.plots(results.filtered)
    return(deconvolution.results.plots(results.filtered))
  })
  names(lgg.permute) <- filter.types.vector
  return(lgg.permute)
}
