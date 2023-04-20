library(dplyr)
library(ggplot2)

set.seed(0)

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
z.data <- c(rnbinom(50, 100, mu = 15), rnbinom(50, 100, mu = 1))
z.data <- c(z.data, rnbinom(50, 100, mu = 1), rnbinom(50, 100, mu = 15))
z.pb <- matrix(z.data, byrow = F, ncol = 2)
s <- s.list.square[[length(s.list.square)]]

#-----------------
# set error params
#-----------------
mean.error = 0
sd.error = 3
z <- z.pb
num.types <- ncol(z)
num.markers <- nrow(z)

#----------------------
# random errors on rows
#----------------------
experiment.label <- "missmatch;row_errors"

# set iterations
iterations.vector <- seq(1000)

# get z series
list.z.error <- lapply(iterations.vector, function(index){
  z.error.rows <- apply(z, 1, function(ri){
    ri + rnorm(1, mean.error, sd.error)
  }) %>% t()
  z.error.rows[z.error.rows < 0] <- 0
  z.error.rows
})

# get simulation results
list.results <- lapply(iterations.vector, function(index){
  simulation_results(s.list.square, z = list.z.error[[index]], 
                     z.pb = z.pb, scale.s = FALSE, 
                     experiment.type = 
                       paste0(experiment.label, ";unscaled"))
})
results.all <- do.call(rbind, list.results)


results.all <- simulation_results(s.list.square, z = z.error.rows, 
                                      z.pb = z.pb, scale.s = FALSE, 
                               experiment.type = paste0(experiment.label, ";unscaled"))

# plot results
ggplot(results.all, aes(x = corr.r.type1, y = abs.error.type1, 
                        color = index, shape = experiment.type)) +
  geom_point() + geom_abline(slope = 1, intercept = 0) + 
  ggtitle("All experiments, composite")

ggplot(results.all, aes(x = p.true.type1, y = p.pred.type1, 
                        shape = experiment.type,
                        color = index)) +
  geom_point() + geom_abline(slope = 1, intercept = 0) + 
  ggtitle("All experiments, composite") + 
  facet_wrap(~experiment.type)

#-------------------------
# random errors on columns
#-------------------------
mean.error <- 10
sd.error <- 10

# set iterations
iterations.vector <- seq(1000)

# get z series
list.z.error <- lapply(iterations.vector, function(index){
  z.error.cols <- apply(z, 2, function(ri){
    ri + rnorm(1, mean.error, sd.error)
  })
  z.error.cols[z.error.cols < 0] <- 0
  z.error.cols
})

# get simulation results
list.results <- lapply(iterations.vector, function(index){
  message(index)
  simulation_results(s.list.square, z = list.z.error[[index]], 
                     z.pb = z.pb, scale.s = FALSE, 
                     experiment.type = 
                       paste0(experiment.label, ";unscaled"))
})
results.all <- do.call(rbind, list.results)


# plot results
ggplot(results.all, aes(x = corr.r.type1, y = abs.error.type1, 
                        color = index, shape = experiment.type)) +
  geom_point() + geom_smooth(method = "glm") + 
  geom_abline(slope = 1, intercept = 0) + 
  ggtitle("All experiments, composite")

ggplot(results.all, aes(x = p.true.type1, y = p.pred.type1, 
                        shape = experiment.type,
                        color = index)) +
  geom_point() + geom_abline(slope = 1, intercept = 0) + 
  ggtitle("All experiments, composite")

#-------------------------
# random errors on entries
#-------------------------
experiment.label <- "missmatch;entries_errors"

# get z series
list.z.error <- lapply(iterations.vector, function(index){
  z.error.entries <- z + rnorm(num.markers*num.types, 
                               mean.error, sd.error)
  z.error.entries[z.error.entries < 0] <- 0
  z.error.entries
})

# get simulation results
list.results <- lapply(iterations.vector, function(index){
  message(index)
  simulation_results(s.list.square, z = list.z.error[[index]], 
                     z.pb = z.pb, scale.s = FALSE, 
                     experiment.type = 
                       paste0(experiment.label, ";unscaled"))
})
results.all <- do.call(rbind, list.results)

# plot results
ggplot(results.all, aes(x = corr.r.type1, y = abs.error.type1, 
                        color = index, shape = experiment.type)) +
  geom_point() + geom_smooth(method = "glm") + 
  geom_abline(slope = 1, intercept = 0) + 
  ggtitle("All experiments, composite")

ggplot(results.all, aes(x = p.true.type1, y = p.pred.type1, 
                        shape = experiment.type,
                        color = index)) +
  geom_point() + geom_abline(slope = 1, intercept = 0) + 
  ggtitle("All experiments, composite")








# get simulation results
results.noscale <- simulation_results(s.list.square, z = z.error.cols, 
                                      z.pb = z.pb, scale.s = FALSE, 
                                      experiment.type = paste0(experiment.label, ";unscaled"))
results.all <- rbind(results.scale, results.noscale)

# plot results
ggplot(results.all, aes(x = corr.r.type1, y = abs.error.type1, 
                        color = experiment.type, shape = experiment.type)) +
  geom_point() + geom_line() + geom_abline(slope = 1, intercept = 0) + 
  ggtitle(experiment.label)

ggplot(results.all, aes(x = p.true.type1, y = p.pred.type1, 
                        color = experiment.type, shape = experiment.type)) +
  geom_point() + geom_abline(slope = 1, intercept = 0) + 
  ggtitle(experiment.label)
