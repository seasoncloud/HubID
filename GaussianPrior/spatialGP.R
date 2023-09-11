# spatial hdp
library(MASS)
library(ggplot2)
#locations 
n = 500 # datapoints

# sampling x and y
sx = runif(n)
sy = runif(n)

s = cbind(sy,sx)
dist = sapply(1:n,function(x) sapply(1:n,function(y) sqrt(sum((s[y,]-s[x,])^2))))

## simulating some data
# possible y

y = rep(1,n)

idx = which(sx < 0.75 & sx > 0.25 & sy >0.25 & sy < 0.75)


y[idx[1:100]] = 3

# putting a little normal noise on.
yobs = y + rnorm(n,0,0.2)
#yobs = y

#weighting 
range = 0.1
phi = 3/range

sigma = exp(-phi*dist)

sigmanorm = diag(1/rowSums(sigma))%*%sigma
weight = c(sigmanorm%*%yobs)

# gaussian mean
range = 1
phi = 3/range

sigma = exp(-phi*dist)

zval = mvrnorm(n=2,yobs,sigma)
gp_mean = colMeans(zval)
gp_median = apply(zval,2,median)

# make it correlated
range = 0.05
phi = 3/range
sigma = exp(-phi*dist)
eig = eigen(sigma)

L = eig$vectors%*%diag(sqrt(eig$values))%*%t(eig$vectors)

ycorr = L%*%yobs

# plotting the results

plotdata = data.frame(s, id = 1:n,yobs,gp_mean,weight,ycorr)
datalong = reshape(plotdata, varying = colnames(plotdata)[-c(1:3)], direction = "long", v.names = 'obs', times = colnames(plotdata)[-c(1:3)])

ggplot(datalong, aes(x = sx, y = sy, col = obs))+
  geom_point(cex = 4, alpha = 0.5)+
  facet_grid(cols = vars(time))





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
