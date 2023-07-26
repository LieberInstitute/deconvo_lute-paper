
#
# This script shows how we use the covariance relationships before/after 
# transformations to random variables.
#


# manually testing covariances
set.seed(10)
x1 <- rnbinom(n = 10, size = 10, mu = 10)  #sample(10)
y1 <- rnbinom(n = 10, size = 10, mu = 10) #sample(10)
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

#--------------------------------------
# transform x2, y2 and compare var, cov
#--------------------------------------
x2 <- x1*10
y2 <- y1*3

# greater covariance suggests greater similarity in variance direction
# (e.g. both variables increase together)
cov(x1, y1) > cov(x2, y2) # TRUE

# covariance magnitude may not tell us strength of relationship, so this may not be helpful
# may consult the correlation coefficient to learn more about this
abs(cov(x1, y1)) > abs(cov(x2, y2)) # FALSE

# strength of relationship increased after transformation (by a small amount)
cor(x1, y1) < cor(x2, y2) # TRUE
cor(x1,y1) - cor(x2, y2) # -5.551115e-17
cor(x1,y1, method = "spearman") - cor(x2, y2, method = "spearman") # 0

var(x2) > var(x1) # TRUE
var(y2) > var(y1) # TRUE












