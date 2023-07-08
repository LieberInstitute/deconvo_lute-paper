source("deconvo_method-paper/code/10_donor-bias-simulations-continued/00_parameters.R")
sapply(libv, library, character.only = T)

results.tables.list <- lapply(dnv, function(dni){
  base.pathi <- file.path(save.dpath, dni, "data")
  fni <- list.files(base.pathi)
  fni <- fni[grepl("^results_table_.", fni)]
  path <- file.path(base.pathi, fni)
  read.csv(path)
})

# make plots
# get plot data
experiment.results.table1 <- results.tables.list[[1]]
experiment.results.table2 <- results.tables.list[[2]]
experiment.results.table1$experiment <- across.batch.label
experiment.results.table2$experiment <- within.batch.label
experiment.results.table <- rbind(experiment.results.table1, 
                                  experiment.results.table2)

# make tall plot data
typev <- unlist(strsplit(experiment.results.table$type_labels[1], ";"))
num.types <- seq(length(typev))
dfp <- do.call(rbind, lapply(num.types, function(ii){
  dfi <- data.frame(method = experiment.results.table$deconvolution_method,
                    experiment = experiment.results.table$experiment,
                    value = experiment.results.table[,paste0("bias.type", ii)])
  dfi$type <- typev[ii]; dfi
}))

# make new plot objects
ggjitter <- ggplot(dfp, aes(x = method, y = value)) + theme_bw() +
  geom_hline(yintercept = 0, color = "red", size = 2) +
  geom_jitter(alpha = 0.5) + geom_boxplot(alpha = 0, color = "cyan", size = 1)
ggjitter + facet_wrap(~experiment+type)









