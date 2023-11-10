
source("model/NMFbatch.R")

# load the data count
# count = read.csv('data/LH_counts.csv', header = T, row.names = 1)
# count = as.matrix(count)

# location = read.csv('data/LH_location.csv', row.names = 1)
# location = as.matrix(location)

path1 = '/gladstone/engelhardt/pelka-collaboration/HuColonCa-FFPE-ImmuOnco-LH_VMSC02001_20220427/HuColonCa-FFPE-ImmuOnco-LH_VMSC02001_20220427/cellpose_cell_by_gene_goodgenes_noblanks.csv'
count = read.csv(path1, row.names = 1)
count = as.matrix(count)
path2 = '/gladstone/engelhardt/pelka-collaboration/HuColonCa-FFPE-ImmuOnco-LH_VMSC02001_20220427/HuColonCa-FFPE-ImmuOnco-LH_VMSC02001_20220427/cellpose_cell_metadata.csv'
metadata = read.csv(path2, row.names = 1)

location = metadata[,c('center_x','center_y')]
location = as.matrix(location)



idx = rowSums(count) > 0

location = location[idx,]
count = count[idx,]

n = nrow(location)

batch_id = groupondist(location,150)

out = nmfspatial_batch(count, 20, location = location, batch = batch_id, lengthscale = 100, initial = 1, smallIter = 500, maxiter = 5000, tolerance = 1e-10)
save(out, file = paste0("modelssaved/cellpose_f20_b150_l100_1e10.RData"))



