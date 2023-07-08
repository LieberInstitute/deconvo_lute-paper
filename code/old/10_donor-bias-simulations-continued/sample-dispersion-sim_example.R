set.seed(0)

num.cells <- 10
genes.per.bin <- 10
m.vector <- seq(1, 30, 1)
c <- 100

mexpr <- t(sapply(m.vector, function(m){
  gamma = c*m/m^2
  rnbinom(genes.per.bin, size = gamma, mu = m)
}))

dfp <- data.frame(mean = rowMeans(mexpr), 
                  var = rowVars(mexpr))
dfp$c <- c

plot(log10(dfp$mean), log10(dfp$var))
abline(0, 1)

library(ggplot2)
c.vector <- seq(0, 1000, 100)
get_dfp <- function(c, num.cells = 10, genes.per.bin = 10, 
                    m.vector = seq(1, 30, 1)){
  mexpr <- t(sapply(m.vector, function(m){
    gamma = c*m/m^2
    rnbinom(genes.per.bin, size = gamma, mu = m)
  }))
  dfp <- data.frame(mean = rowMeans(mexpr), var = rowVars(mexpr))
  dfp$c <- c
  return(dfp)
}

dfp.new <- do.call(rbind, lapply(c.vector, get_dfp))
dfp.new$c <- as.character(dfp.new$c)

ggplot(dfp.new, aes(x = mean, y = var, group = c, color = c)) +
  geom_smooth(se = F) + geom_abline(slope = 1, intercept = 0) +
  scale_x_log10() + scale_y_log10()

ggplot(dfp.new, aes(x = mean, y = var, group = c, color = c)) +
  geom_smooth() + geom_abline(slope = 1, intercept = 0) +
  geom_point(alpha = 0.5)

ggplot(dfp.new, aes(x = c, y = var, group = c, color = c)) +
  geom_violin(draw_quantiles = 0.5)

