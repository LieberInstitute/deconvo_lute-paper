
#
# This script shows how we use the covariance relationships before/after 
# transformations to random variables.
#

libv <- c("ggplot2")
sapply(libv, library, character.only = T)
set.seed(10)

#-----------------------
# get starting variables
#-----------------------
x1 <- rnbinom(n = 10, size = 10, mu = 10)  #sample(10)
y1 <- rnbinom(n = 10, size = 10, mu = 10) #sample(10)

# terms for manually testing covariances
#mean.x1 <- mean(x1)
#mean.y1 <- mean(y1)
#N <- length(x1)

#------------------------------
# get covariance of 2 variables
#------------------------------
cov(x1,y1) # -3.7

#-----------------------------
# get variance of one variable
#-----------------------------
cov(x1, x1) == var(x1) 
var(x1, x1) == var(x1)

#---------------------------------------------------------------
# transform x2, y2 and compare var, cov -- static transformation
#---------------------------------------------------------------
x2 <- x1*10
y2 <- y1*3

# greater covariance suggests greater VARIANCE DIRECTION SIMILARITY (not relation strength)
# (e.g. both variables increase together)
cov(x1, y1) < 0 # TRUE
cov(x2, y2) < 0 # TRUE
cov(x1, y1) > cov(x2, y2) # TRUE

# covariance magnitude does not tell us strength of relationship, 
# rather is just the geometric mean of variances in common for x,y
abs(cov(x1, y1)) > abs(cov(x2, y2)) # FALSE

# strength of relationship increased after transformation (by a small amount)
cor(x1, y1) < cor(x2, y2) # TRUE
cor(x1,y1) - cor(x2, y2) # -5.551115e-17
cor(x1,y1, method = "spearman") - cor(x2, y2, method = "spearman") # 0

var(x2) > var(x1) # TRUE
var(y2) > var(y1) # TRUE

#----------------------------------------------------------------------------
# transform x1, y1 and compare var, cov, cor -- define dynamic transformation
#----------------------------------------------------------------------------
# can we transform x and y such that their correlation decreases, and their covariance remains negative?
#

# For Y = 3, X = 10
# 
# we have:
# x1/y1 = 1/1 = 1
# X/Y = 10/3 = 3.3333...
# 10/3 = x2/y2 = (x1*X)/(y1*Y)
# 10/3 = x2/y2
# 
# solve for X
# x2/y2 = (x1*X)/(y1*Y)
# X = (x2*y1*Y)/(y2*x1) = fract * y1*Y/(y2*x1)
# 
# solve for Y
# x2/y2 = (x1*X)/(y1*Y)
# (y1*Y)(x2/y2) = x1*X
# Y = (y2*x1*X)/(x2*y1) = 1/fract * (x1*X/y1)
# 
# Note: if x1 == y1, expect: X.adj/Y.adj == fract == X/Y
#

adj.weight.x <- function(WeightX, WeightY, x1, y1){
  # x1 : first cell type
  # y1 : second cell type
  fract <- WeightX/WeightY
  fract*y1*WeightY/x1 %>% round()
}
get.adj <- function(x1, y1, WeightX, WeightY){
  # pseudoalgorithm
  # given a vector of normalizations, get the fraction, x/y
  # now adjust vectors x1, y1 such that the ratio is achieved
  X.adj <- adj.weight.x(WeightX = X, WeightY = Y, x1 = x1, y1 = y1)
  Y.adj <- adj.weight.x(WeightX = Y, WeightY = X, x1 = x1, y1 = y1)
  list(x.adj=x1*X.adj, y.adj=y1*Y.adj)
}

# for identical terms (x1==y1)
X.adj <- adj.weight.x(WeightX = X, WeightY = Y, x1 = x1[1], y1 = y1[1])
Y.adj <- adj.weight.x(WeightX = Y, WeightY = X, x1 = x1[1], y1 = y1[1])
X.adj # 10
Y.adj # 3

# inspect static relation
# starting normalization factors
fract <- 10/3
fract # 3.333
# get static relation -- when terms identical
iter.index <- 1
x1[iter.index] # 9
y1[iter.index] # 9
x2[iter.index] # 90
y2[iter.index] # 27
x2[iter.index]/y2[iter.index] # 3.333
# get static relation -- when terms nonidentical
iter.index <- 2
x1[iter.index] # 8
y1[iter.index] # 7
x2[iter.index] # 80
y2[iter.index] # 21
x2[iter.index]/y2[iter.index] # 3.809524
# summarize fract across x, y
summary(x1/y1)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.5000  0.7500  0.8889  1.2309  1.1321  3.7500 
summary(x2/y2)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.667   2.500   2.963   4.103   3.774  12.500
# plot sets
plot(x1, y1)
plot(x2, y2)
plot(x1-x2, y1-y2)
ggplot(data.frame(x1=x1,y1=y1),aes(x=x1,y=y1)) + geom_point() + geom_abline()
ggplot(data.frame(x2=x2,y2=y2),aes(x=x2,y=y2)) + geom_point() + geom_abline() + xlim(0,max(x2)) + ylim(0,max(y2))
ggplot(data.frame(x1m2=(x1-x2),y1m2=(y1-y2)),aes(x=x1m2,y=y1m2)) + geom_point() + geom_abline(slope=1,intercept=0)

# inspect dynamic relation
adj.xy <- get.adj(x1, y1, WeightX, WeightY)
x2 <- adj.xy[[1]]
y2 <- adj.xy[[2]]
X <- 10
Y <- 3
fract <- X/Y
summary(x2/y2) # summary of the fractions x/y
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.667   2.500   2.963   4.103   3.774  12.500
# comparison to original fract
round(x2/y2,4)==round(fract,4) # round for 3.333...~1/3
# [1]  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# plot relationship
plot(x1, y1)
plot(x2, y2)
plot(x1-x2, y1-y2)
dfp <- rbind(data.frame(x = x1, y = y1, type = rep("unadj", length(x1))),
             data.frame(x = x2, y = y2, type = rep("adj", length(x2))))
ggplot(dfp, aes(x = x, y = y, color = type)) + geom_point() + facet_wrap(~type) +
  geom_abline(slope = 1, intercept = 0)

#--------------------------
# compare static vs dynamic
#--------------------------
dfp <- dfp
dfp[dfp$type=="adj",]$type <- "dynamic"
dfp <- rbind(dfp, data.frame(x = x1*X, y = y1*Y, type = rep("static", length(x1))))
ggplot(dfp, aes(x = x, y = y, color = type)) + geom_point() + facet_wrap(~type) +
  geom_abline(slope = 1, intercept = 0)

# cor tests
# unadj
cor(dfp[dfp$type=="unadj",]$x,dfp[dfp$type=="unadj",]$y) # -0.2988827
# static
cor(dfp[dfp$type=="static",]$x,dfp[dfp$type=="static",]$y) # -0.2988827
# dynamic
cor(dfp[dfp$type=="dynamic",]$x,dfp[dfp$type=="dynamic",]$y) # [1] 0.9354588

# cov tests
# unadj
cov(dfp[dfp$type=="unadj",]$x,dfp[dfp$type=="unadj",]$y) # -3.7
# static
cov(dfp[dfp$type=="static",]$x,dfp[dfp$type=="static",]$y) # -111
# dynamic
cov(dfp[dfp$type=="dynamic",]$x,dfp[dfp$type=="dynamic",]$y) # 1093.662



