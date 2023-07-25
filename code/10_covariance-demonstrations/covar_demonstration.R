# manually testing covariances
set.seed(10)
x <- sample(10)
y <- sample(10)
mean.x <- mean(x)
mean.y <- mean(y)
N <- length(x)

# attempt1
term.vector <- ((x-mean.x)*(y-mean.y))/N
sum(term.vector)

# attempt2
sum((x*y)-(mean.x*mean.y))

# attempt3
cov(x,y)