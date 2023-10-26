library(Rcpp)
library(RcppArmadillo)
sourceCpp("model/NMFspatial2.cpp")

path = "/gladstone/engelhardt/pelka-collaboration/HuColonCa-FFPE-ImmuOnco-LH_VMSC02001_20220427/HuColonCa-FFPE-ImmuOnco-LH_VMSC02001_20220427"
setwd(path)

data = read.csv('cellpose_cell_by_gene_goodgenes_noblanks.csv', row.names = 1)
data = as.matrix(data)


out = nmfgen(data, 10, iter = 100)
save(out, file = paste0("modelssaved/cellpose_result.RData"))


