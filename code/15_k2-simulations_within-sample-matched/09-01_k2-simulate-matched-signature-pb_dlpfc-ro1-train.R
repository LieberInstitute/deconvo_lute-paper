
# setup simulation
# perform deconvolution
# uses matched pseudobulk and signature matrices within samples
# includes null cell sizes adjustment

source("deconvo_method-paper/code/15_k2-simulations_within-sample-matched/00_parameters-script-set-15.R")
sapply(libv, library, character.only = T)

# load marker expression data
se.path <- here("deconvo_method-paper/outputs/15_k2-simulations_within-sample-matched/",
                        "se-markers-overlap-filters_dlpfc-ro1-train.rda")
se <- get(load(se.path))

# filter overlapping markers
markers.data <- metadata(se)$marker.data$markers.overlaps
markers.high.overlap <- lapply(markers.data, function(data){
  data[data[,2] >= 5, 1]
}) %>% unlist()
# get overall mean ratios
markers.se <- lute(sce = SingleCellExperiment(se), 
                   celltype.variable = "k2", 
     deconvolution.algorithm = NULL,
     markers.per.type = 100)

se <- se[rownames(se) %in% markers.high.overlap,]
dim(se)
length(intersect(rownames(se), markers.data$neuron[,1])) # 49
length(intersect(rownames(se), markers.data$glial[,1])) # 49

# simulations -- group-matched z and y-pseudobulk
# experiment parameters
assay.name <- "counts"
celltype.variable <- "k2"
sample.id.variable <- "Sample"
S.pb <- c("glial" = 3, "neuron" = 10)
deconvolution.algorithm <- "nnls"
experiment <- deconvolution.experiment(y = NULL, sce = se, s = S.pb, 
                                       assay.name = assay.name,
                                       sample.id.variable = sample.id.variable, 
                                       celltype.variable = celltype.variable, 
                                       markers.vector = NULL)
# save
experiment.name <- "experiment_within-slide-matched_dlpfc-ro1-train.rda"
experiment.path <- here("deconvo_method-paper", "outputs", experiment.name)
save(experiment, file = experiment.path)
# plots
experiment$plots.list$neuron$proportions.scatterplot
experiment$plots.list$neuron$abs.error.barplot
experiment$plots.list$neuron$abs.error.jitterbox

# simulations -- missmatched z, y-pseudobulk
# experiment parameters
# get true, pred proportions for sample pairs
z.series <- c(unique(se[["Sample"]]), "all")
sample.id.vector <- unique(se[["Sample"]])
list.results.all <- lapply(z.series, function(z.id){
  message(z.id)
  z.set <- se
  if(!z.id == "all"){z.set <- se[,se[["Sample"]]==z.id]}
  z <- signature_matrix_from_sce(z.set, "k2")
  amount.true.z <- table(z.set$k2)
  names(amount.true.z) <- paste0("amount.true.z.", names(amount.true.z))
  prop.true.z <- amount.true.z %>% prop.table()
  names(prop.true.z) <- paste0("prop.true.z.", names(prop.true.z))
  group.y.vector <- sample.id.vector[!sample.id.vector == z.id]
  result.list <- lapply(group.y.vector, function(sample.id.ypb){
    # validate a group result
    y <- ypb_from_sce(se[,se[["Sample"]]==sample.id.ypb], "counts", "k2", S = S.pb)
    # with adjustment
    result.adj <- nnlsParam(z = z, y = y, s = S.pb) %>% deconvolution()
    prop.pred.adj <- result.adj/sum(result.adj)
    result.unadj <- nnlsParam(z = z, y = y, s = c("glial" = 1, "neuron" = 1)) %>% 
      deconvolution()
    prop.pred.unadj <- result.unadj/sum(result.unadj)
    # prop true
    amount.true.y <- table(se[,se[["Sample"]]==sample.id.ypb]$k2)
    prop.true.y <- amount.true.y %>% prop.table()
    # append names
    names(amount.true.y) <- paste0("amount.true.y.", names(amount.true.y))
    names(prop.pred.adj) <- paste0("prop.pred.adj.", names(prop.pred.adj))
    names(prop.pred.unadj) <- paste0("prop.pred.unadj.", names(prop.pred.unadj))
    names(prop.true.z) <- paste0("prop.true.z.", names(prop.true.z))
    names(prop.true.y) <- paste0("prop.true.y.", names(prop.true.y))
    c(z.id, sample.id.ypb, amount.true.z, amount.true.y, 
      prop.true.z, prop.true.y, prop.pred.adj, prop.pred.unadj)
  })
  result.table <- do.call(rbind, result.list) %>% as.data.frame()
  colnames(result.table)[1:2] <- c("sample.id.z", "sample.id.y")
  for(c in c(3:14)){result.table[,c] <- as.numeric(result.table[,c])}
  result.table
})
results.all.table <- do.call(rbind, list.results.all) %>% as.data.frame()
for(c in c(3:14)){results.all.table[,c] <- as.numeric(results.all.table[,c])}
results.all.table$error.adj.neuron <- results.all.table$prop.pred.adj.neuron - 
  results.all.table$prop.true.y.neuron
results.all.table$abs.error.adj.neuron <- abs(results.all.table$error.adj.neuron)
results.all.table$error.unadj.neuron <- results.all.table$prop.pred.unadj.neuron - 
  results.all.table$prop.true.y.neuron
results.all.table$abs.error.unadj.neuron <- abs(results.all.table$error.unadj.neuron)
#results.all.table$difference.true.prop.neuron <- 
#  results.all.table$prop.true.z.neuron-results.all.table$prop.true.y.neuron

# plots -- errors
# error values
ggplot(results.all.table, aes(x = error.unadj.neuron, y = error.adj.neuron)) +
  geom_abline(slope = 1, intercept = 0, color = "red") + 
  geom_vline(xintercept = 0, alpha = 0.5) +
  geom_hline(yintercept = 0, alpha = 0.5) +
  geom_point() + facet_wrap(~sample.id.z)
# absolute errors
ggplot(results.all.table, 
       aes(x = abs.error.unadj.neuron, y = abs.error.adj.neuron)) +
  geom_abline(slope = 1, intercept = 0, color = "red") + 
  geom_vline(xintercept = 0, alpha = 0.5) +
  geom_hline(yintercept = 0, alpha = 0.5) +
  geom_point() + facet_wrap(~sample.id.z) + 
  xlim(0, 0.8) + ylim(0, 0.8)
# plots -- proportions
plot.table <- do.call(rbind, 
                      lapply(c("prop.pred.adj.neuron",
                               "prop.pred.unadj.neuron"),
                             function(column.name){
                                      data.frame(pred = results.all.table[,column.name],
                                                 type = column.name,
                                                 true = 
                                                   results.all.table[,"prop.true.y.neuron"])
}))
ggplot(plot.table, aes(x = true, y = pred)) +
  geom_point() + geom_abline(slope = 1, intercept = 0) + facet_wrap(~type)

# plots -- amount differences by errors
ggplot(plot.table, aes(x = abs.difference, y = abs.error)) +
  geom_point() + geom_abline(slope = 1, intercept = 0) + facet_wrap(~type) +
  facet_wrap(~z.id)

# plots -- amount z neuron by errors
ggplot(results.all.table, aes(x = amount.true.z.neuron, y = abs.error.adj.neuron)) +
  geom_point() + geom_abline(slope = 1, intercept = 0) + facet_wrap(~type) +
  facet_wrap(~sample.id.z) + scale_x_log10()
