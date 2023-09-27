source("./scripts/09_fast/00_param.R")
mae <- get(load("./outputs/01_mae/mae_analysis_append.rda"))
mae <- mae[,colData(mae)$sample.id=="Br8492_mid",]
dfs.iter <- dfs.series(s.glial.series = seq(1, 40, 2))
df.res <- multigroup_bias_matched(sample.id.vector = "Br8492_mid", 
                        df.true.list = metadata(mae[[1]])[["list.df.true.k2"]], 
                        y.unadj = mae[["bulk.rnaseq"]][,1], 
                        dfs = dfs.iter, 
                        sce = mae[[1]], 
                        y.group.name = "batch.id2",
                        celltype.variable = "k2", 
                        assay.name = "logcounts")
df.res <- dfres_postprocess(df.res)
summary(df.res$error.neuron)
min(df.res$error.neuron)
library(ggplot2)
ggplot(df.res, aes(x = "sample", y = df.res$error.neuron)) + geom_violin(draw_quantiles=0.5)
filter.error <- df.res$error.neuron==min(df.res$error.neuron)
s.optimal <- c(df.res[filter.error,]$s.glial,df.res[filter.error,]$s.neuron)
names(s.optimal) <- c("glial", "neuron")
s.optimal
source("./source/00_dataset_summaries.R")
source("./source/00_deconvo_plots.R")
list.plots.res <- get_dfres_plots(df.res)
list.plots.res$heatmaps$heatmap1
list.plots.res$heatmaps$heatmap6