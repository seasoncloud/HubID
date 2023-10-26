# Function for batch updates of spatial nmf

nmfspatial_batch = function(data, noSignatures, location, batch, maxiter = 10000, tolerance = 1e-8, initial = 5, smallIter = 100){
    weights = list()
    batch_list = list()
    for(i in unique(batch)){
        index = which(batch == i)
        batch_list[[i]] = index

        X = location[index,]
        X2 = rowSums(X^2)
        X2 = matrix(X2, nrow = length(X2), ncol = length(X2), byrow = F)
        r = X2 - 2*X%*%t(X) + t(X2)
        if(any(r<0)){
            warning("Some distances were smaller than zero! Try scaling up the locations.")
            r[r<0] = 0
        }
        
        dist = sqrt(r)

        weight[[i]] = dist
    }

    out = list
    out$b = batch_list
    out$w = weights

    return(out)

    
}