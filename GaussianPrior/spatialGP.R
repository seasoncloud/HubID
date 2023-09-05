# spatial hdp
library(MASS)
library(ggplot2)
#locations 
n = 100 # datapoints
#sx = c(runif(n/8,0,0.05),runif(n/8,0.6,1),runif(n/2,0.6,1),runif(n/4,0,0.4))
#sy = c(runif(n/8,0,0.05),runif(n/8,0.6,1),runif(n/2,0,0.4),runif(n/4,0.6,1))

# sampling x and y
sx = runif(n)
sy = runif(n)

s = cbind(sy,sx)

# creating correlation matrix
phi = 0.002

#range
3/phi

sigma = sapply(1:n,function(x) sapply(1:n,function(y) exp(-phi*sqrt(sum((s[y,]-s[x,])^2)))))

Tdata = 1000
#cdf_mu = rbeta(1,1,2)

#mu = qnorm(cdf_mu)
mu = 0
zval = mvrnorm(n=1,rep(mu,n),sigma)
x = rnorm(100)
plotdata = data.frame(col = zval,s, id = 1:n, x = x)

ggplot(plotdata, aes(x=sx,y = sy, col = zval, label = id))+
  geom_point(cex = 2)
  geom_text(hjust = 0, vjust = 0)

# possible y
y = rep(1,n)

idx = which(sx < 0.75 & sx > 0.25 & sy >0.25 & sy < 0.75)


y[idx[1:12]] = 1.1

zval_spatial = exp(colMeans(mvrnorm(n=100,y,sigma)))
zval_independent = exp(colMeans(mvrnorm(n=100,y,diag(n))))


plotdata = data.frame(zval,zval_ind,s, id = 1:n)

ggplot(plotdata, aes(x=sx,y = sy, col = zval, label = id))+
  geom_point(cex = 2)

ggplot(plotdata, aes(x=sx,y = sy, col = zval_ind, label = id))+
  geom_point(cex = 2)
