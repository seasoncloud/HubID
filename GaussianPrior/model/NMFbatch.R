# Function for batch updates of spatial nmf
library(Rcpp)
library(RcppArmadillo)
sourceCpp("model/NMFspatial2.cpp")

nmfspatial_batch = function(data, noSignatures, location, lengthscale, batch, maxiter = 10000, tolerance = 1e-8, initial = 5, smallIter = 100){
    weights = list()
    batch_list = list()
    for(i in unique(batch)){
        index = which(batch == i)
        batch_list[[i]] = index - 1

        X = location[index,]
        X2 = rowSums(X^2)
        X2 = matrix(X2, nrow = length(X2), ncol = length(X2), byrow = F)
        r = X2 - 2*X%*%t(X) + t(X2)
        if(any(r<0)){
            warning("Some distances were smaller than zero! Try scaling up the locations.")
            r[r<0] = 0
        }
        
        dist = sqrt(r)
        
        # calculating covariance
        sigma = exp(-dist/lengthscale)
        sigma[sigma < 0.5] = 0

        weights[[i]] = sigma/rowSums(sigma)
    }

    out = nmfspatialbatch(data = data, noSignatures = noSignatures, weight = weights, batch = batch_list, maxiter = maxiter, tolerance = tolerance, initial = initial, smallIter = smallIter)

    return(out)  
}