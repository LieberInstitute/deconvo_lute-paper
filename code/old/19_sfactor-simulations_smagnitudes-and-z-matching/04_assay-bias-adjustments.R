library(dplyr)
library(ggplot2)
set.seed(0)

#-----------------
# helper functions
#-----------------
simulation_results <- function(s.list.experiment, z, 
                               y.adjustment.type = c("pb", "z", NULL),
                               z.pb = NULL,
                               scale.s = FALSE, p.true = c(0.2, 0.8),
                               experiment.type = "s_squared_relation"){
  if(is(z.pb, "NULL")){z.pb <- z}
  corr.type1 <- corr.type2 <- "NA"
  if(nrow(z)>2){
    corr.type1 <- cor.test(z[,1], z.pb[,1], method = "pearson")$estimate
    corr.type2 <- cor.test(z[,2], z.pb[,2], method = "pearson")$estimate  
  }
  y.series <- lapply(s.list.experiment, function(s){
    zs <- sweep(z.pb, 2, s, "*")
    ypb <- t(t(p.true) %*% t(zs))
    ypb
  })
  index.vector <- seq(length(y.series))
  experiment.list <- lapply(index.vector, function(index){
    y.experiment <- y.series[[index]]
    s.iteration <- s.list.experiment[[index]]
    z.iteration <- z
    if(scale.s){z.iteration <- sweep(z.iteration, 2, s.iteration, "*")}
    if(is(y.adjustment.type, "NULL")){
      message("proceeding without y adjustment")
    } else if(y.adjustment.type == "pb"){
      ypb.z <- t(t(p.true) %*% t(z.iteration))
      lm.data <- cbind(y.experiment, ypb.z) %>% as.data.frame()
      fit.z <- lm(y.experiment ~ ypb.z, data = lm.data)
      y.experiment <- fit.z$fitted.values
    } else if(y.adjustment.type == "z"){
      lm.data <- cbind(y.experiment, z.iteration) %>% as.data.frame()
      z.variables <- paste0("z", seq(num.types))
      colnames(lm.data) <- c("y.experiment", z.variables)
      fit.text <- paste0("lm(y.experiment ~ ", 
                         paste0(z.variables, collapse = " + "),
                         ", data = lm.data)")
      fit.z <- eval(parse(text = fit.text))
      y.experiment <- fit.z$fitted.values
    } else{
      stop("error, unidentified adjustment type.")
    }
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

adjustment_simulations <- function(z, s.list, z.pb = NULL, p.true.list = list(c(0.2, 0.8))){
  list.simulations <- lapply(p.true.list, function(p.true.iter){
    # without scaling
    result.unscale1 <- simulation_results(s.list.experiment = s.list, 
                                          z = z, 
                                          z.pb = z.pb, 
                                          y.adjustment.type = "pb", 
                                          scale.s = FALSE,
                                          experiment.type = "adj_pb;unscaled",
                                          p.true = p.true.iter)
    result.unscale2 <- simulation_results(s.list.experiment = s.list, 
                                          z = z, 
                                          z.pb = z.pb, 
                                          y.adjustment.type = "z", 
                                          scale.s = FALSE,
                                          experiment.type = "adj_z;unscaled",
                                          p.true = p.true.iter)
    result.unscale3 <- simulation_results(s.list.experiment = s.list, 
                                          z = z, 
                                          z.pb = z.pb, 
                                          y.adjustment.type = NULL, 
                                          scale.s = FALSE,
                                          experiment.type = "unadj;unscaled",
                                          p.true = p.true.iter)
    # with scaling
    result.scale1 <- simulation_results(s.list, z, z.pb = zs, "pb", scale.s = TRUE, 
                                        experiment.type = "adj_pb;scale",
                                        p.true = p.true.iter)
    result.scale2 <- simulation_results(s.list, z, z.pb = zs, "z", scale.s = TRUE,
                                        experiment.type = "adj_z;scale",
                                        p.true = p.true.iter)
    result.scale3 <- simulation_results(s.list, z, z.pb = zs, NULL, scale.s = TRUE,
                                        experiment.type = "unadj;scale",
                                        p.true = p.true.iter)
    # get plot data
    result.scale <- rbind(result.scale1, result.scale2, result.scale3)
    result.unscale <- rbind(result.unscale1, result.unscale2, result.unscale3)
    result.all <- rbind(result.scale, result.unscale)
    return(result.all)
  })
  sim.results <- do.call(rbind, list.simulations) %>% as.data.frame()
  plot.data <- data.frame(prop.pred.type1 = sim.results$p.pred.type1,
                          prop.true.type1 = sim.results$p.true.type1,
                          abs.error.type1 = sim.results$abs.error.type1,
                          corr.r.type1 = sim.results$corr.r.type1,
                          experiment.type = sim.results$experiment.type,
                          s.type1 = sim.results$s.type1)
  plot1 <- ggplot(plot.data, aes(x = prop.true.type1, y = prop.pred.type1,
                                 color = s.type1)) + geom_point() +
    geom_abline(slope = 1, intercept = 0) + facet_wrap(~experiment.type) +
    xlim(0, 0.3) + ylim(0, 0.3)
  plot2 <- ggplot(plot.data, aes(x = corr.r.type1, y = abs.error.type1,
                                 color = s.type1)) + 
    geom_point() + facet_wrap(~experiment.type)
  plot3 <- ggplot(plot.data, aes(x = experiment.type, y = abs.error.type1)) + 
    geom_jitter() + geom_boxplot(alpha = 0, color = "cyan")
  lgg <- list(plot.data = plot.data, plot1 = plot1, plot2 = plot2, plot3 = plot3)
  list.results <- list(results = sim.results, lgg = lgg)
  return(list.results)
}

make_z <- function(mu.low = 1, mu.high = 15, seed.num = 0){
  set.seed(seed.num)
  z.data <- c(rnbinom(50, 100, mu = mu.high), 
              rnbinom(50, 100, mu = mu.low))
  z.data <- c(z.data, rnbinom(50, 100, mu = mu.low), 
              rnbinom(50, 100, mu = mu.high))
  z <- matrix(z.data, byrow = F, ncol = 2)
}

#--------------------------
# get main simulation terms
#--------------------------
mean.error = 0
sd.error = 3
num.types <- 2
num.markers <- 100

# get s
s <- c(3, 10)
# get p.true
p.true <- c(0.2, 0.8)

# get z
z.data <- c(rnbinom(50, 100, mu = 15), rnbinom(50, 100, mu = 1))
z.data <- c(z.data, rnbinom(50, 100, mu = 1), rnbinom(50, 100, mu = 15))
z <- matrix(z.data, byrow = F, ncol = 2)
# get pb
y.pb1 <- t(t(p.true) %*% t(z))
zs <- sweep(z, 2, s, "*")
y.pb2 <- t(t(p.true) %*% t(zs))

#----------------------------------------------------------
# test adjustment, scaling for a series of true proportions
#----------------------------------------------------------
s.vector <- seq(1, 10, 1)
s.list.square <- lapply(s.vector, function(item){c(item, item^2)})

# for single p.true set
results <- adjustment_simulations(z = z, z.pb = zs, s.list = s.list.square)

# for multiple p.true sets
p.true.type1.vector <- seq(0.05, 0.35, 0.01)
p.true.list <- lapply(p.true.type1.vector, function(value){c(value,1-value)})
results <- adjustment_simulations(z = z, z.pb = zs, s.list = s.list.square, 
                                  p.true.list = p.true.list)

#-----------------------------------------------------------------
# test adjustment, scaling for a series of randomized z references
#-----------------------------------------------------------------

# get missmatched z based on correlation
mu.low <- 1
mu.high <- 15
iter.vector <- seq(1e6)
z.new.start <- make_z(iter.vector[1])
corr.r <- cor.test(z.new.start[,1], z[,1], method = "spearman")$estimate
while(corr.r > 0.5){
  message(iter.new)
  iter.new <- iter + 1
  mu.low = mu.low * 0.99
  mu.high = mu.high * 1.01
  z.new.start <- make_z(mu.low = mu.low, mu.high = mu.high,
                        seed.num = iter.vector[iter.new])
  corr.r <- cor.test(z.new.start[,1], z[,1], method = "spearman")$estimate
}

corr.r.start <- make_z(iter[1])




z.series <- lapply(seq(100), function(index){make_z(index)})
list.results <- lapply(z.series, function(z.iteration){
  adjustment_simulations(z = z.iteration, s.list = list(c(3, 10)))
})

plot.data <- do.call(rbind, lapply(list.results, function(results){results$lgg$plot.data}))
  
plot1 <- ggplot(plot.data, aes(x = prop.true.type1, y = prop.pred.type1,
                               color = s.type1)) + geom_point() +
  geom_abline(slope = 1, intercept = 0) + facet_wrap(~experiment.type) +
  xlim(0, 0.3) + ylim(0, 0.3)

plot2 <- ggplot(plot.data, aes(x = corr.r.type1, y = abs.error.type1,
                               color = s.type1)) + 
  geom_point() + geom_smooth(method = 'lm') + facet_wrap(~experiment.type)

plot3 <- ggplot(plot.data, aes(x = experiment.type, y = abs.error.type1)) + 
  geom_jitter() + geom_boxplot(alpha = 0, color = "cyan")  


#--------------------------------------------------------------------------------------
# test adjustment, scaling for a series of true proportions and randomized z references
#--------------------------------------------------------------------------------------



#-------------------------------------
# adjustment -- using only pseudobulks
#-------------------------------------
lm.data <- cbind(y.pb2, y.pb1) %>% as.data.frame()
colnames(lm.data) <- c("y.pb.zs", "y.pb.z")
fit.pb <- lm(y.pb.zs ~ y.pb.z, data = lm.data)
# plot results
plot(density(fit.pb$residuals))
plot.data <- data.frame(original = lm.data$y.pb.zs, 
                        fitted = fit.pb$fitted.values)
ggplot(plot.data, aes(x = original, y = fitted)) + 
  geom_point() + geom_abline(slope = 1, intercept = 0)
# corr compare
cor.test(y.pb1, fit.pb$fitted.values, method = "spearman")$estimate
cor.test(y.pb1, y.pb2, method = "spearman")$estimate

#-------------------------------------
# adjustment -- using signature matrix
#-------------------------------------
lm.data <- cbind(y.pb2, z) %>% as.data.frame()
colnames(lm.data) <- c("ypb.zs", "z1", "z2")
fit.z <- lm(ypb.zs ~ z1 + z2, data = lm.data)
# plot results
plot(density(fit.z$residuals))
plot.data <- data.frame(original = lm.data$ypb.zs, fitted = fit$fitted.values)
ggplot(plot.data, aes(x = original, y = fitted)) + 
  geom_point() + geom_abline(slope = 1, intercept = 0)
# corr compare
cor.test(y.pb1, fit.z$fitted.values, method = "spearman")$estimate
cor.test(y.pb1, y.pb2, method = "spearman")$estimate

#---------------------------------------
# deconvolution before and after fitting
#---------------------------------------
list.y.data <- list(
  unadj = y.pb2,
  adj.pb = fit.pb$fitted.values,
  adj.z = fit.z$fitted.values
)

list.result <- lapply(list.y.data, function(yi){
  pred.iter <- nnls::nnls(z, yi)$x
  pred.iter/sum(pred.iter)
})
result <- do.call(rbind, list.result) %>% as.data.frame()
colnames(result) <- c("pred.type1", "pred.type2")
result$y.type <- rownames(result)
result$abs.error.type1 <- abs(result$pred.type1 - p.true[1])

#-----------------------------------------------------------
# sfactor improvement with and without pseudobulk adjustment
#-----------------------------------------------------------
p.true <- c(0.2, 0.8)
s <- c(3, 10)
# get z
z.data <- c(rnbinom(50, 100, mu = 15), rnbinom(50, 100, mu = 1))
z.data <- c(z.data, rnbinom(50, 100, mu = 1), rnbinom(50, 100, mu = 15))
z <- z <- matrix(z.data, byrow = F, ncol = 2)

# get pb
y.pb1 <- t(t(p.true) %*% t(z))
zs <- sweep(z.pb, 2, s, "*")
y.pb2 <- t(t(p.true) %*% t(zs))




