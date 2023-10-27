library(Rcpp)
library(RcppArmadillo)
source("model/NMFbatch.R")

# load the data count
count = read.csv('data/LH_counts.csv', header = T, row.names = 1)
count = as.matrix(count)

location = read.csv('data/LH_location.csv', row.names = 1)
location = as.matrix(location)


out = nmfspatial_batch(count, 10,location = location, batch = rep(1,nrow(location)), lengthscale = 1000, initial = 1, smallIter = 100, maxiter = 200)
save(out, file = paste0("modelssaved/LH_f10_trybatch.RData"))



