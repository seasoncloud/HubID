library(Rcpp)
library(RcppArmadillo)
sourceCpp("model/NMFspatial2.cpp")

# load the data count
count = read.csv('data/LH_counts.csv', header = T, row.names = 1)
count = as.matrix(count)

location = read.csv('data/LH_location.csv', row.names = 1)
location = as.matrix(location)

X = location
X2 = rowSums(X^2)
X2 = matrix(X2, nrow = length(X2), ncol = length(X2), byrow = F)
dist = sqrt(X2 - 2*X%*%t(X) + t(X2))


sigma = exp(-0.005*dist)
sigma[sigma < 0.5] = 0
sigma = sigma/rowSums(sigma)



out = nmfspatial(count, 10,weight = sigma, initial = 1, smallIter = 100, maxiter = 5000)
save(out, file = paste0("modelssaved/LH5K_f10_l0001.RData"))



