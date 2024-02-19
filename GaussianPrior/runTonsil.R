library(Rcpp)
library(RcppArmadillo)
source("model/NMFbatch.R")

# load the data count
count = read.csv('data/data_mf/tonsil_protein.csv', header = T, row.names = 1)
count = as.matrix(count)

range(rowSums(count))

location = read.csv('data/data_mf/tonsil_location.csv', row.names = 1)
location = as.matrix(location)

r1 = range(location[,1])
r2 = range(location[,2])
density = (r1[2] - r1[1])*(r2[2] - r2[1])/nrow(count)
density = round(density,0)
density = density/10
print(density)
density = 20
batch_id = groupondist(location,size = 20000)

# start = Sys.time()
# out = nmfspatial_batch(count, 6, location = location, batch = batch_id, lengthscale = density, initial = 1, maxiter = 1000, tolerance = 1e-8, error_freq = 10, kernel_cutoff = 0.2, normalize = T)
# stop = Sys.time()
# out$time = stop - start
# print(stop - start)

# save(out, file = paste0("modelssaved/tonsil_f6_l",density,"_norm.RData"))


row_sum = rowSums(count)
count = count/row_sum*mean(row_sum) 

start = Sys.time()
out = nmfgen(count, 6, iter = 1000)
stop = Sys.time()
out$time = stop - start
print(stop - start)
save(out, file = paste0("modelssaved/tonsil_f6_nmfgen_norm_i1000.RData"))


