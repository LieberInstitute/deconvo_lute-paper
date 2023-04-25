library(dplyr)
library(ggplot2)

set.seed(0)

#-----------------
# helper functions
#-----------------
make_experiment_list <- function(z, s = c(3, 10), p = c(0.2, 0.8), z.pb = NULL){
  num.experiment.s <- num.experiment.z <- 
    num.experiment.p <- num.experiment.zpb <- 1
  if(is(s, "list")){
    num.experiment.s <- length(s)
  }
  if(is(z, "list")){
    num.experiment.z <- length(z)
  }
  if(is(p, "list")){
    num.experiment.p <- length(p)
  }
  if(is(z.pb, "list")){
    num.experiment.zpb <- length(z.pb)
  }
  num.experiment <- max(
    c(num.experiment.s, num.experiment.z, num.experiment.p, num.experiment.zpb))
  lapply(seq(num.experiment), function(iter){
    if(is(s, "list")){
      s.iter <- s[[iter]]
    } else{
      s.iter <- s
    }
    if(is(z, "list")){
      z.iter <- z[[iter]]
    } else{
      z.iter <- z
    }
    if(is(p, "list")){
      p.iter <- p[[iter]]
    } else{
      p.iter <- p
    }
    if(is(z.pb, "list")){
      zpb.iter <- z.pb[[iter]]
    } else{
      zpb.iter <- z.pb
    }
    zpb.s.iter <- sweep(zpb.iter, 2, s.iter, "*")
    y.iter <- t(t(p.iter) %*% t(zpb.s.iter))
    
    # get z correlations
    seq.types <- seq(ncol(z.iter))
    type.z.corr <- lapply(seq.types, function(type.iter){
      cor.test(z.iter[,type.iter], zpb.iter[,type.iter],
               method = "pearson")$estimate
    }) %>% unlist()
    names(type.z.corr) <- paste0("type", seq.types)
    
    list(index = iter, s = s.iter, z = z.iter, 
         p = p.iter, z.pb = z.pb, y = y.iter,
         z.corr = type.z.corr)
  })
}

experiment_results <- function(list.experiment.info){
  index.vector <- length(list.experiment.info) %>% seq()
  experiment.list <- lapply(index.vector, function(index){
    experiment.info <- list.experiment.info[[index]]
    y.iter <- experiment.info[["y"]]
    s.iter <- experiment.info[["s"]]
    p.iter <- experiment.info[["p"]]
    z.corr <- experiment.info[["z.corr"]]
    z.iter <- experiment.info[["z"]]
    z.iter.scale <- sweep(z.iter, 2, s.iter, "*")
    
    # unscaled results
    p.pred.unscale <- nnls::nnls(A = z.iter, b = y.iter)$x
    p.pred.unscale <- p.pred.unscale/sum(p.pred.unscale)
    
    # scaled results
    p.pred.scale <- nnls::nnls(A = z.iter.scale, b = y.iter)$x
    p.pred.scale <- p.pred.scale/sum(p.pred.scale)
    
    # errors
    error.unscale.vector <- p - p.pred.unscale
    error.scale.vector <- p - p.pred.scale
    
    # new row data
    new.row <- c(index, p.iter, s.iter, z.corr, p.pred.unscale, p.pred.scale,
                 error.unscale.vector, error.scale.vector,
                 abs(error.unscale.vector), abs(error.scale.vector))
    seq.types <- length(s.iter) %>% seq()
    names(new.row) <- c("index", 
                        paste0("p.true.type", seq.types),
                        paste0("s.type", seq.types),
                        paste0("z.corr.type", seq.types),
                        paste0("p.pred.unscaled.type", seq.types),
                        paste0("p.pred.scaled.type", seq.types),
                        paste0("error.unscaled.type", seq.types),
                        paste0("error.scaled.type", seq.types),
                        paste0("abs.error.unscaled.type", seq.types),
                        paste0("abs.error.scaled.type", seq.types))
    new.row
  })
  experiment.table <- do.call(rbind, experiment.list) %>% as.data.frame()
  # colnames(experiment.table) <- c("index", "s.type1", "s.type2", "p.pred.type1", "p.pred.type2")
  return(experiment.table)
}

simulation_results <- function(z, s = c(3, 10), p = c(0.2, 0.8), z.pb = NULL,
                               experiment.type = "new_experiment"){
  require(dplyr)
  # get experiment list info
  if(is(z.pb, "NULL")){
    if(is(z, "list")){
      z.pb <- z[[1]]
    } else{
      z.pb <- z
    }
  }
  list.experiment.info <- make_experiment_list(s = s, z = z, p = p, z.pb = z.pb)
  # get results table
  experiment.results.table <- experiment_results(list.experiment.info)
  experiment.results.table$experiment.type <- experiment.type
  return(experiment.results.table)
}

#---------------------
# set z with transform
#---------------------
# set error params
mean.error = 0
sd.error = 3
z <- z.pb
num.types <- ncol(z)
num.markers <- nrow(z)


# set z
z.data <- c(rnbinom(50, 100, mu = 15), rnbinom(50, 100, mu = 1))
z.data <- c(z.data, rnbinom(50, 100, mu = 1), rnbinom(50, 100, mu = 15))
z.pb <- matrix(z.data, byrow = F, ncol = 2)

# random errors on rows
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

#-----------------------
# get simulation results
#-----------------------
# get series across varied p values
p1 <- seq(0.1, 0.4, 1e-2)
list.results <- lapply(p1, function(p1.iter){
  p2.iter <- 1-p1.iter
  simulation_results(z = list.z.error, p = c(p1.iter, p2.iter))
})
results.all <- do.call(rbind, list.results)

# make tall data for plotting
results.tall1 <- data.frame(p.true.type1 = results.all$p.true.type1,
                            p.pred.type1 = results.all$p.pred.scaled.type1,
                            z.corr.type1 = results.all$z.corr.type1,
                            abs.error.type1 = results.all$abs.error.scaled.type1)
results.tall1$type <- "scaled"
results.tall2 <- data.frame(p.true.type1 = results.all$p.true.type1,
                            p.pred.type1 = results.all$p.pred.unscaled.type1,
                            z.corr.type1 = results.all$z.corr.type1,
                            abs.error.type1 = results.all$abs.error.unscaled.type1)
results.tall2$type <- "unscaled"
results.tall <- rbind(results.tall1, results.tall2)

#------------------
# make scatterplots
#------------------
# scatterplots -- p.true vs. p.predicted, by scale
ggplot(results.tall, aes(x = p.true.type1, y = p.pred.type1, color = z.corr.type1)) +
  geom_point(alpha = 0.4) + geom_abline(slope = 1, intercept = 0) + 
  ggtitle("All experiments, composite") + geom_smooth(method = "glm") +
  facet_wrap(~type)

# scatterplots -- z.corr vs. abs.error, by scale
ggplot(results.tall, aes(x = z.corr.type1, y = abs.error.type1)) +
  geom_point(alpha = 0.2) + geom_abline(slope = 1, intercept = 0) + 
  ggtitle("All experiments, composite") + geom_smooth(method = "glm") +
  facet_wrap(~type)
