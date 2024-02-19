# Function for batch updates of spatial nmf
library(Rcpp)
library(RcppArmadillo)
sourceCpp("model/NMFspatial2.cpp")

dist_fun = function(X){
    X2 = rowSums(X^2)
    X2 = matrix(X2, nrow = length(X2), ncol = length(X2), byrow = F)
    r = X2 - 2*X%*%t(X) + t(X2)
    if(any(r<0)){
        warning("Some distances were smaller than zero! Try scaling up the locations.")
        r[r<0] = 0
    }
        
    return(r)
}

dist_index = function(X,index){
    X2 = rowSums(X^2)
    r = sum(X[index,]^2) - 2*X[index,]%*%t(X) + X2
    if(any(r<0)){
        warning("Some distances were smaller than zero! Try scaling up the locations.")
        r[r<0] = 0
    } 
    return(r)
}

groupondist = function(location, size = NULL, no_groups = NULL){
    n = nrow(location)
    left = c(1:n)

    if(is.null(no_groups) & is.null(size)){
        stop("You must determine the size or number of groups")
    }

    if(is.null(size)){
        size = ceiling(n/no_groups)
    }

    i = 1
    batch_vec = rep(paste0("b",0),n)
    while(length(left) > size){ 
        start = sample(length(left),1)
        dist = dist_index(location[left,],start)
        batch_index = order(dist)[1:size]
        batch_vec[left[batch_index]] = paste0("b",i)
        i = i+1
        
        left = left[-batch_index]
    }
    return(batch_vec)
}

nmfspatial_batch = function(data, noSignatures, location, lengthscale, batch = 1, maxiter = 10000, tolerance = 1e-8, initial = 5, smallIter = 100, error_freq = 10,kernel_cutoff = 0.5,normalize = TRUE){
    
    if(normalize){
        row_sum = rowSums(data)
        data = data/row_sum*mean(row_sum)
    }
    
    unique_batches = unique(batch)
    if(length(unique_batches) == 1){
        print("Everything is run in one batch")

        dist = dist_fun(location)
        
        # calculating covariance
        sigma = exp(-dist/(lengthscale^2))
        sigma[sigma < kernel_cutoff] = 0

        weight = sigma/rowSums(sigma)

        out = nmfspatial(data = data, noSignatures = noSignatures, weight = weight, maxiter = maxiter, tolerance = tolerance, initial = initial, smallIter = smallIter)

    }else{
        weights = list()
        batch_list = list()
        for(i in unique_batches){
            index = which(batch == i)
            batch_list[[i]] = index - 1

            X = location[index,]
        
            dist = dist_fun(X)
        
            # calculating covariance
            sigma = exp(-dist/(lengthscale^2))
            sigma[sigma < kernel_cutoff] = 0

            weights[[i]] = sigma/rowSums(sigma)
        }

        if(initial == 1){
            out = nmfspatialbatch2(data = data, noSignatures = noSignatures, weight = weights, batch = batch_list, maxiter = maxiter, tolerance = tolerance, error_freq = error_freq)
        }else{
            out = nmfspatialbatch(data = data, noSignatures = noSignatures, weight = weights, batch = batch_list, maxiter = maxiter, tolerance = tolerance, initial = initial, smallIter = smallIter, error_freq = error_freq)
        }

    }
    
    return(out)
}