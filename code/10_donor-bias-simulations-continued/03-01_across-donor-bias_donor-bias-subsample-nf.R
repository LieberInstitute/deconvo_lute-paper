#!/usr/bin/env R

# Author: Sean Maden
#
# Run donor bias subsampling experiments.
#
#

source("deconvo_method-paper/code/10_donor-bias-simulations-continued/00_parameters.R")
sapply(libv, library, character.only = T)

# set up a new subsample experiment
sce <- get(load(sce.list.path))[[celltype.variable]]; rm(lscef)
proportions <- sce[[k.marker.variable]] %>% table() %>% prop.table()
lsub <- prepare_subsample_experiment(sce, scale.factor = scale.factor,
                                     iterations = iterations, groups.per.iteration = 3,
                                     method.vector = methodv, celltype.variable = celltype.variable,
                                     group.variable = group.variable, assay.name = assay.name,
                                     cell.proportions = proportions, count.minimum = count.min,
                                     seed.num = seed.num, which.save = c("sce", "tp", "ypb", "li"),
                                     save.fnstem = save.fnstem, base.directory.workflow.data = "data", verbose = TRUE)
setwd(workflow.path)
# save workflow table
workflow.table <- lsub$workflow.table
write.csv(workflow.table, file = workflow.table.path)
# save table iterations
index.vector <- seq(1, nrow(workflow.table), num.batch)
workflow.table.filename.ter <- paste0("workflow-table_",save.fnstem,".csv")
for(iteration.number in seq(length(index.vector))){
  index.start <- index.vector[iteration.number]
  index.end <- (index.start+num.batch-1)
  if(index.end > nrow(workflow.table)){index.end <- nrow(workflow.table)}
  workflow.tablei <- workflow.table[index.start:index.end,]
  workflow.tablei.fname.iter <- workflow.table.fnamei <- paste0("workflow-table-iter", 
                                        iter.num, "_", save.fnstem,".csv")
  workflow.tablei.fpath <- file.path(base.directory.workflow.data, workflow.tablei.fname.iter)
  write.csv(workflow.tablei, file = workflow.tablei.fpath, row.names = F)
  message("finished with iter ", iter.num)
}

# run the workflow
# navigate to:
# /deconvo_method-paper/outputs/10_donor-bias-simulations-continued/r-nf_deconvolution_donor-bias
#
# run:
# bash ./sh/r-nf_run_workflow.table-iter.sh

# analyze results table
results.filt <- "results_table_.*"
data.dpath <- file.path(save.dpath, rnf.dname, "data")
lfv <- list.files(data.dpath)
rt.fname <- lfv[grepl(results.filt, lfv)]
rt.fpath <- file.path(save.dpath, rnf.dname, "data", rt.fname)
rt <- read.csv(rt.fpath)

# prep results table
varv <- c("prop.pred.type1", "prop.pred.type2", "bias.type1", "bias.type2", "rmse.types")
method.varname <- "deconvolution_method"
# filter duplicate entries
rtf <- rt
rtf$label <- paste0(rtf$iterations_index, ";", rtf$deconvolution_method)
rtf <- rtf[!duplicated(rtf$label),]
methodv <- unique(rtf[,method.varname])

# jitter plots
lgg <- lapply(varv, function(varname){
  dfp <- rtf; dfp$value <- dfp[,varname]
  ggplot(dfp, aes(x = deconvolution_method, y = value)) +
    geom_jitter(alpha = 0.5) + ggtitle(varname) +
    geom_boxplot(alpha = 0, color = "cyan", size = 1)
})
# bias plots
plot1 <- lgg[[3]] + theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  geom_hline(yintercept = 0, alpha = 1, color = "red")
plot2 <- lgg[[4]] + theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0, alpha = 1, color = "red")
jpeg("ggjitterbox-comp-bias_inter-sample-bias.jpg", 
     width = 6, height = 5, units = "in", res = 400)
grid.arrange(plot1, plot2, nrow = 2,
             layout_matrix = matrix(c(1,1,1,2,2,2,2), ncol = 1),
             bottom = "Deconvolution method")
dev.off()

# plot pairs
# varname <- "prop.pred.type1"
# iterate on plot variables
lgg <- lapply(varv, function(varname){
  dfp <- do.call(cbind, lapply(methodv, function(methodi){
    filt <- rtf[,method.varname]==methodi; rti <- rtf[filt,]
    rti <- rti[order(rti$iterations_index),]; rti[,varname]
  }))
  dfp <- as.data.frame(dfp); colnames(dfp) <- methodv
  ggpairs(dfp) + ggtitle(varname)
})
names(lgg) <- varv

lgg$prop.pred.type1
lgg$prop.pred.type2
lgg$bias.type1
lgg$bias.type2
lgg$rmse.types

# violin plots
lgg <- lapply(varv, function(varname){
  dfp <- rtf; dfp$value <- dfp[,varname]
  ggplot(dfp, aes(x = deconvolution_method, y = value)) +
    geom_violin(draw_quantiles = 0.5) + ggtitle(varname)
})
jpeg("ggvp-comp_inter-sample-bias.jpg", width = 5, height = 4,
     units = "in", res = 400)
grid.arrange(lgg[[1]], lgg[[2]], lgg[[3]], lgg[[4]], lgg[[5]],
             layout_matrix = matrix(c(1,2,3,4,5), nrow = 2))
dev.off()
# bias plots
plot1 <- lgg[[3]] + theme_bw() +
  theme(axis.text.x = element_blank(),
                          axis.title.x = element_blank()) +
  geom_hline(yintercept = 0, alpha = 0.2)
plot2 <- lgg[[4]] + theme_bw() +
  theme(axis.title.x = element_blank(),
                          axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0, alpha = 0.2)
jpeg("ggvp-comp-bias_inter-sample-bias.jpg", 
     width = 6, height = 5, units = "in", res = 400)
grid.arrange(plot1, plot2, nrow = 2,
             layout_matrix = matrix(c(1,1,1,2,2,2,2), ncol = 1),
             bottom = "Deconvolution method")
dev.off()

# scatter plot bias
# no facet
ggpt <- ggplot(rtf, aes(x = bias.type1, y = bias.type2, 
                group = deconvolution_method, 
                color = deconvolution_method)) +
  geom_point(alpha = 0.2) + geom_abline(slope = 1, intercept = 0) +
  theme_bw() + 
  geom_hline(yintercept = 0, alpha = 0.5, color = "gray") +
  geom_vline(xintercept = 0, alpha = 0.5, color = "gray")
# facet
ggpt <- ggplot(rtf, aes(x = bias.type1, y = bias.type2)) +
  geom_point(alpha = 0.2) + geom_abline(slope = 1, intercept = 0) +
  theme_bw() + 
  geom_hline(yintercept = 0, alpha = 0.5, color = "gray") +
  geom_vline(xintercept = 0, alpha = 0.5, color = "gray")
ggpt <- ggpt + facet_wrap(~deconvolution_method)
jpeg('ggpt-facet-method_inter-sample-bias.jpg', 
     width = 6, height = 6, 
     units = "in", res = 400)
ggpt
dev.off()
