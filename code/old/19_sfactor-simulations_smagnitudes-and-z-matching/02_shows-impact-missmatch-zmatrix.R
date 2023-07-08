library(dplyr)
library(ggplot2)

#-----------------
# helper functions
#-----------------
simulation_results <- function(s.list.experiment, z, z.pb = NULL,
                               scale.s = FALSE, p.true = c(0.2, 0.8),
                               experiment.type = "s_squared_relation"){
  if(is(z.pb, "NULL")){z.pb <- z}
  corr.type1 <- corr.type2 <- "NA"
  if(nrow(z)>2){
    corr.type1 <- cor.test(z[,1], z.pb[,1], method = "pearson")$estimate
    corr.type2 <- cor.test(z[,2], z.pb[,2], method = "pearson")$estimate  
  }
  y.series <- lapply(s.list.experiment, function(s){
    zs <- sweep(z.pb, 2, s, "*"); t(t(p.true) %*% t(zs))
  })
  index.vector <- seq(length(y.series))
  experiment.list <- lapply(index.vector, function(index){
    y.experiment <- y.series[[index]]
    s.iteration <- s.list.experiment[[index]]
    z.iteration <- z
    if(scale.s){z.iteration <- sweep(z.iteration, 2, s.iteration, "*")}
    p.predicted <- nnls::nnls(A = z.iteration, b = y.experiment)$x
    p.predicted <- p.predicted/sum(p.predicted)
    c(index, s.iteration, p.predicted)
  })
  experiment.table <- do.call(rbind, experiment.list) %>% as.data.frame()
  colnames(experiment.table) <- c("index", "s.type1", "s.type2", "p.pred.type1", "p.pred.type2")
  experiment.table$corr.r.type1 <- round(corr.type1, 3)
  experiment.table$corr.r.type2 <- round(corr.type2, 3)
  experiment.table$s.difference <- experiment.table$s.type1 - experiment.table$s.type2
  experiment.table$s.abs.difference <- abs(experiment.table$s.difference)
  experiment.table$p.true.type1 <- p.true[1]
  experiment.table$p.true.type2 <- p.true[2]
  experiment.table$error.type1 <- experiment.table$p.pred.type1 - experiment.table$p.true.type1
  experiment.table$error.type2 <- experiment.table$p.pred.type2 - experiment.table$p.true.type2
  experiment.table$abs.error.type1 <- abs(experiment.table$error.type1)
  experiment.table$abs.error.type2 <- abs(experiment.table$error.type2)
  experiment.table$experiment.type <- experiment.type
  return(experiment.table)
}

#-----------
# get s list
#-----------
s.vector <- seq(1, 10, 1)
s.list.square <- lapply(s.vector, function(item){c(item, item^2)})

#---------------------
# set z with transform
#---------------------
# set z
z.pb <- matrix(c(0, 1, 0, 1, 1, 0), byrow = T, ncol = 2)
s <- s.list.square[[length(s.list.square)]]
z <- sweep(z, 2, s, "*")

# results labels
results1.label <- "transform;matched"
results2.label <- "notransform;matched"
results3.label <- "transform;missmatch"
results4.label <- "notransform;missmatch"

# matched, with transform
results1 <- simulation_results(s.list.square, z = z, 
                               scale.s = TRUE, experiment.type = results1.label)
# matched, no transform
results2 <- simulation_results(s.list.square, z = z, 
                               scale.s = FALSE, experiment.type = results2.label)
# missmatched, with transform
results3 <- simulation_results(s.list.square, z = z, z.pb = z.pb, scale.s = TRUE, 
                               experiment.type = results3.label)
# missmatched, no transform
results4 <- simulation_results(s.list.square, z = z, z.pb = z.pb, scale.s = FALSE, 
                               experiment.type = results4.label)

results.all <- rbind(results1, results2, results3, results4)

# plots
ggplot(results1, aes(x = s.abs.difference, y = abs.error.type1)) +
  geom_point() + geom_line() + geom_abline(slope = 1, intercept = 0) +
  ggtitle(results1.label)

ggplot(results2, aes(x = s.abs.difference, y = abs.error.type1)) +
  geom_point() + geom_line() + geom_abline(slope = 1, intercept = 0) +
  ggtitle(results2.label)

ggplot(results3, aes(x = s.abs.difference, y = abs.error.type1)) +
  geom_point() + geom_line() + geom_abline(slope = 1, intercept = 0) +
  ggtitle(results3.label)

ggplot(results4, aes(x = s.abs.difference, y = abs.error.type1)) +
  geom_point() + geom_line() + geom_abline(slope = 1, intercept = 0) +
  ggtitle(results4.label)

ggplot(results.all, aes(x = s.abs.difference, y = abs.error.type1)) +
  geom_point() + geom_line() + geom_abline(slope = 1, intercept = 0) + 
  facet_wrap(~experiment.type) + ggtitle("All experiments, faceted")

ggplot(results.all, aes(x = p.true.type1, y = p.pred.type1, 
                        color = experiment.type, shape = experiment.type)) +
  geom_point() + geom_line() + geom_abline(slope = 1, intercept = 0) + 
  ggtitle("All experiments, composite") + facet_wrap(~experiment.type)

ggplot(results.all, aes(x = corr.r.type1, y = abs.error.type1, 
                        color = experiment.type, shape = experiment.type)) +
  geom_point() + geom_line() + geom_abline(slope = 1, intercept = 0) + 
  ggtitle("All experiments, composite") + facet_wrap(~experiment.type)
