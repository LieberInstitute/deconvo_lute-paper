#----------
# load data
#----------

df1.fname <- "results-table_inter-sample_100iter-present.csv"
df1 <- read.csv(df1.fname)
df2.fname <- "results-table_intra-sample_100iter-present.csv"
df2 <- read.csv(df2.fname)

#-------------
# bias by type
#-------------
dfp1 <- data.frame(bias = c(df1$bias.type1, df1$bias.type2),
                   type = c(rep("glial", nrow(df1)), rep("neuron", nrow(df1))))
dfp1$expt <- "across"
dfp1$method <- df1$deconvolution_method
dfp2 <- data.frame(bias = c(df2$bias.type1, df2$bias.type2),
                   type = c(rep("glial", nrow(df2)), rep("neuron", nrow(df2))))
dfp2$expt <- "within"
dfp2$method <- df1$deconvolution_method
dfp <- rbind(dfp1, dfp2)

ggplot(dfp, aes(x = expt, y = bias)) + ylab("Bias") + xlab("Experiment") +
  facet_wrap(~type+method) + theme_bw() +
  geom_hline(yintercept = 0, color = "blue") + geom_jitter(alpha = 1e-1) +
  geom_violin(draw_quantiles = 0.5, color = "black", alpha = 0)

ggplot(dfp, aes(x = expt, y = bias)) + ylab("Bias") + xlab("Experiment") +
  facet_wrap(~type+method) + theme_bw() +
  geom_hline(yintercept = 0, color = "black") + geom_jitter(alpha = 1e-1) +
  geom_boxplot(color = "blue", alpha = 0) + 
  ggtitle("Bias (true - predicted proportions)")

#--------------------------------
# plot -- scatter plot by method
#--------------------------------
df1$experiment <- "across"
df2$experiment <- "within"
dfr <- rbind(df1, df2)

# scatter plot -- bias
metric.plot <- title.str <- 'bias'
# type predictions by method
dfp1 <- dfr[dfr$deconvolution_method=="nnls",]
dfp2 <- dfr[dfr$deconvolution_method=="music",]
dfp1$iter.id <- paste0(dfp1$iterations_index, ";", dfp1$experiment)
dfp2$iter.id <- paste0(dfp2$iterations_index, ";", dfp2$experiment)
row1.orderv <- order(match(dfp1$iter.id, dfp2$iter.id))
dfp1 <- dfp1[row1.orderv,]
cond <- identical(dfp1$iterations_index, dfp2$iterations_index)
# get plot data
typev <- c(".type1", ".type2")
names(typev) <- unique(unlist(strsplit(dfp1$type_labels, ";")))
dfp <- do.call(rbind, lapply(c(".type1", ".type2"), function(typei){
  var.str <- paste0(metric.plot, typei)
  dfpi <- data.frame(nnls = dfp1[,var.str], music = dfp2[,var.str])
  dfpi$experiment <- dfp1$iter.id
  dfpi$type <- names(typev[typev==typei]); dfpi
}))

dfp$id <- gsub(".*;", "", dfp$experiment)
ggplot(dfp, aes(x = nnls, y = music)) + 
  geom_abline(slope = 1, intercept = 0, color = "black") +
  geom_hline(yintercept = 0, color = "gray") + geom_vline(xintercept = 0, color = "gray") +
  geom_point(alpha = 0.5) +
  facet_wrap(~type+id) + theme_bw() + geom_smooth() +
  ggtitle("Bias (true - predicted proportions)")

