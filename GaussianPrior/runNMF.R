library(Rcpp)
library(RcppArmadillo)
sourceCpp("model/NMFspatial2.cpp")

path1 = '/gladstone/engelhardt/pelka-collaboration/HuColonCa-FFPE-ImmuOnco-LH_VMSC02001_20220427/HuColonCa-FFPE-ImmuOnco-LH_VMSC02001_20220427/cellpose_cell_by_gene_goodgenes_noblanks.csv'
count = read.csv(path1, row.names = 1)
count = as.matrix(count)



# load the data count
#data = read.csv('data/LH_counts.csv', header = T, row.names = 1)
#data = as.matrix(data)
idx = rowSums(count > 0) >= 20


count = count[idx,]

out = nmfgen(count, 25, iter = 3000)

save(out, file = paste0("modelssaved/cellpose_f25_nmfgen_over20genes.RData"))


