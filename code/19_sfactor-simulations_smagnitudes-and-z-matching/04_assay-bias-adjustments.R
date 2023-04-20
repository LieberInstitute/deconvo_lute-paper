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

#-----------
# get s list
#-----------
s.vector <- seq(1, 10, 1)
s.list.square <- lapply(s.vector, function(item){c(item, item^2)})

# without scaling
result.unscale1 <- simulation_results(s.list.square, z, "pb", scale.s = FALSE)
result.unscale2 <- simulation_results(s.list.square, z, "z", scale.s = FALSE)
result.unscale3 <- simulation_results(s.list.square, z, NULL, scale.s = FALSE)

# with scaling
result.scale1 <- simulation_results(s.list.square, z, "pb", scale.s = TRUE)
result.scale2 <- simulation_results(s.list.square, z, "z", scale.s = TRUE)
result.scale3 <- simulation_results(s.list.square, z, NULL, scale.s = TRUE)

#---------------------
# set z with transform
#---------------------
# set z
z.data <- c(rnbinom(50, 100, mu = 15), rnbinom(50, 100, mu = 1))
z.data <- c(z.data, rnbinom(50, 100, mu = 1), rnbinom(50, 100, mu = 15))
z.pb <- matrix(z.data, byrow = F, ncol = 2)
s <- s.list.square[[length(s.list.square)]]

#--------------------------
# get main simulation terms
#--------------------------
mean.error = 0
sd.error = 3
num.types <- ncol(z)
num.markers <- nrow(z)

# get s
s <- c(3, 10)
# get p.true
p.true <- c(0.2, 0.8)

# get z
z.data <- c(rnbinom(50, 100, mu = 15), rnbinom(50, 100, mu = 1))
z.data <- c(z.data, rnbinom(50, 100, mu = 1), rnbinom(50, 100, mu = 15))
z <- z <- matrix(z.data, byrow = F, ncol = 2)
# get pb
y.pb1 <- t(t(p.true) %*% t(z))
zs <- sweep(z.pb, 2, s, "*")
y.pb2 <- t(t(p.true) %*% t(zs))

#-----------------------------
# define adjustment as medians
#-----------------------------
# quantify bias
yz.bind <- cbind(z, y.pb)
row.wise.diffs <- apply(yz.bind, 1, function(ri){
  abs(ri[1:2]-ri[3]) %>% median()
})
bias <- median(row.wise.diffs)
zm <- z+bias

# summary plots
plot(density(row.wise.diffs))
plot.data <- data.frame(z1 = z[,1], zm1 = zm[,1])
ggplot(plot.data, aes(x = z1, y = zm1)) + geom_point() + 
  geom_abline(slope = 1, intercept = 1) + xlim(0, 100) + ylim(0,100)

# evaluate adjustment
corr.r.type1.start <- cor.test(y.pb, z[,1])$estimate
corr.r.type1.end <- cor.test(y.pb, zm[,1])$estimate

corr.r.type2.start <- cor.test(y.pb, z[,1])$estimate
corr.r.type2.end <- cor.test(y.pb, zm[,1])$estimate

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




