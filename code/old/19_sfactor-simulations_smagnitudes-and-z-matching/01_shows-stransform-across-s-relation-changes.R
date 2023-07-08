library(dplyr)
library(ggplot2)

#-----------------
# helper functions
#-----------------
simulation_results <- function(s.list.experiment, z, 
                               scale.s = FALSE, p.true = c(0.2, 0.8),
                               experiment.type = "s_squared_relation"){
  y.series <- lapply(s.list.experiment, function(s){
    zs <- sweep(z, 2, s, "*"); t(t(p.true) %*% t(zs))
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

#------------
# set s lists
#------------
s.vector <- seq(1, 10, 1)
# x2 relationship
s.list.x2 <- lapply(s.vector, function(item){c(item, item + item)})
# squared relationship
s.list.square <- lapply(s.vector, function(item){c(item, item^2)})
# shared magnitude relationship
s.list.magnitude <- lapply(s.vector, function(item){c(item, item + 1)})

#------
# set z
#------
z <- matrix(c(0, 1, 1, 0), nrow = 2)

#-----------------------------
# results without scaling on s
#-----------------------------
results1.notransform <- simulation_results(s.list.x2, 
                                           z = z, scale.s = FALSE, 
                                           experiment.type = "s_x2_relation")
results2.notransform <- simulation_results(s.list.square, 
                                           z = z, scale.s = FALSE, 
                                           experiment.type = "s_squared_relation")
results3.notransform <- simulation_results(s.list.magnitude, 
                                           z = z, scale.s = FALSE, 
                                           experiment.type = "s_magnitude_relation")
results.all.notransform <- rbind(results1.notransform, results2.notransform, 
                                 results3.notransform)

ggplot(results.all.notransform, aes(x = s.abs.difference, y = abs.error.type1)) +
  geom_point() + geom_line() + geom_abline(slope = 1, intercept = 0) +
  facet_wrap(~experiment.type) + ggtitle("No transform")

#--------------------------
# results with scaling on s
#--------------------------
z <- matrix(c(0, 1, 1, 0), nrow = 2)
results1.withtransform <- simulation_results(s.list.x2, 
                                           z = z, scale.s = TRUE, 
                                           experiment.type = "s_x2_relation")
results2.withtransform <- simulation_results(s.list.square, 
                                           z = z, scale.s = TRUE, 
                                           experiment.type = "s_squared_relation")
results3.withtransform <- simulation_results(s.list.magnitude, 
                                           z = z, scale.s = TRUE, 
                                           experiment.type = "s_magnitude_relation")
results.all.withtransform <- rbind(results1.withtransform, 
                                 results2.withtransform, 
                                 results3.withtransform)

ggplot(results.all.withtransform, aes(x = s.abs.difference, y = abs.error.type1)) +
  geom_point() + geom_line() + geom_abline(slope = 1, intercept = 0) +
  facet_wrap(~experiment.type) + ggtitle("With transform")

# plot matched results
plot.data <- data.frame(no.transform = results.all.notransform$abs.error.type1,
                        with.transform = results.all.withtransform$abs.error.type1,
                        experiment.type = results.all.withtransform$experiment.type)

ggplot(plot.data, aes(x = no.transform, y = with.transform)) +
  geom_point() + geom_line() + 
  geom_abline(slope = 1, intercept = 0, col = "red", lty = 2) +
  facet_wrap(~experiment.type) + ggtitle("Type1 abs. error")

data1 <- results.all.notransform
data2 <- results.all.withtransform
data1$transform <- "no"
data2$transform <- "yes"
plot.data <- rbind(data1, data2) %>% as.data.frame()
ggplot(plot.data, aes(x = transform, y = abs.error.type1)) + 
  geom_jitter() + geom_boxplot(col = "cyan", alpha = 0)

ggplot(plot.data, aes(x = transform, y = abs.error.type1)) + 
  geom_jitter() + geom_boxplot(col = "cyan", alpha = 0) + 
  facet_wrap(~experiment.type)

ggplot(plot.data, aes(x = p.true.type1, y = p.pred.type1)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0) + facet_wrap(~experiment.type)
