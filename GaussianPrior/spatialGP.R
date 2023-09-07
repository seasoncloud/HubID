# spatial hdp
library(MASS)
library(ggplot2)
#locations 
n = 1000 # datapoints
#sx = c(runif(n/8,0,0.05),runif(n/8,0.6,1),runif(n/2,0.6,1),runif(n/4,0,0.4))
#sy = c(runif(n/8,0,0.05),runif(n/8,0.6,1),runif(n/2,0,0.4),runif(n/4,0.6,1))

# sampling x and y
sx = runif(n)
sy = runif(n)

s = cbind(sy,sx)
dist = sapply(1:n,function(x) sapply(1:n,function(y) sqrt(sum((s[y,]-s[x,])^2))))

# creating correlation matrix
phi = 30

#range
range = 0.03
phi = 3/range

sigma = exp(-phi*dist)


## simulating some data
# possible y

y = rep(1,n)

idx = which(sx < 0.75 & sx > 0.25 & sy >0.25 & sy < 0.75)


y[idx] = 2

yobs = y + rnorm(n,0,0.5)

sigmanorm = diag(1/rowSums(sigma))%*%sigma
znew = c(sigmanorm%*%yobs)

#zval_spatial = colMeans(mvrnorm(n=10,yobs,sigma))

#zval_independent = colMeans(mvrnorm(n=10,yobs,diag(n)))


plotdata = data.frame(yobs,znew, s, id = 1:n)
datalong = reshape(plotdata, varying = colnames(plotdata)[c(1:2)], direction = "long", v.names = 'obs', times = c('original observation','spatial information through gp'))

ggplot(datalong, aes(x = sx, y = sy, col = obs))+
  geom_point(cex = 2)+
  facet_grid(cols = vars(time))

length(znew)
dim(s)
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



ggplot(plotdata, aes(x=sx,y = sy, col = zval_spatial, label = id))+
  geom_point(cex = 2)

ggplot(plotdata, aes(x=sx,y = sy, col = zval_independent, label = id))+
  geom_point(cex = 2)
rgamma(100,1,1)
