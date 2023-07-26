
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

#---------------------------------------------------------------------
# transform x1, y1 and compare var, cov, cor -- dynamic transformation
#---------------------------------------------------------------------
# can we transform x and y such that their correlation decreases, and their covariance remains negative?
#

# For Y = 3, X = 10
# 
# 10/3 = x2/y2 = (x1*X)/(y1*Y)
# 10/3 = x2/y2
# 
# solve for X
# x2/y2 = (x1*X)/(y1*Y)
# X = (x2*y1)/(y2*Y*x1) = fract * y1/(Y*x1)
# 
# solve for Y
# x2/y2 = (x1*X)/(y1*Y)
# (y1*Y)(x2/y2) = x1*X
# Y = (y2*x1*X)/(x2*y1) = 1/fract * (x1*X/y1)
# 
# Note: if x1 == y1, expect: Y.adj/X.adj == fract == X/Y
#

adj.weight.x <- function(Y, x1, y1, fract){
  fract*y1/(Y*x1) %>% round()
}

adj.weight.y <- function(X, x1, y1, fract){
  (fract^-1)*(x1*X)/y1 %>% round()
}
X.adj <- adj.weight.x(Y = 3, x1 = x1[1], y1 = y1[1], fract = 10/3)
Y.adj <- adj.weight.y(X = 10, x1 = x1[1], y1 = y1[1], fract = 10/3)


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







# pseudoalgorithm
# given a vector of normalizations, get the fraction, x/y
fract <- 10/3
# now adjust vectors x1, y1 such that the ratio is achieved
x.iter1 <- x1[1]
y.iter1 <- y1[1]
x.iter2 <- x.iter1*fract
y.iter2 <- y.iter1/fract
x2 <- x1*10
y2 <- y1*3









