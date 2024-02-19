
source("model/NMFbatch.R")

# load dataset 1
path1 = '/gladstone/engelhardt/pelka-collaboration/HuColonCa-FFPE-ImmuOnco-LH_VMSC02001_20220427/HuColonCa-FFPE-ImmuOnco-LH_VMSC02001_20220427/cellpose_cell_by_gene_goodgenes_noblanks.csv'
count = read.csv(path1, row.names = 1)
count1 = as.matrix(count)
path2 = '/gladstone/engelhardt/pelka-collaboration/HuColonCa-FFPE-ImmuOnco-LH_VMSC02001_20220427/HuColonCa-FFPE-ImmuOnco-LH_VMSC02001_20220427/cellpose_cell_metadata.csv'
metadata = read.csv(path2, row.names = 1)

location = metadata[,c('center_x','center_y')]
location1 = as.matrix(location)

idx = rowSums(count1) >= 20
write.csv(idx,file='data/cellidx_LH.csv',row.names=F)

location1 = location1[idx,]
count1 = count1[idx,]

batch_id1 = groupondist(location1,size = 5000, lengthscale = 100)
batch_id1 = paste0("1_",batch_id1)

# load dataset 2
path1 = '/gladstone/engelhardt/pelka-collaboration/HumanColonCancerPatient1/cell_by_gene.csv'
count = read.csv(path1, row.names = 1)
count2 = as.matrix(count)
path2 = '/gladstone/engelhardt/pelka-collaboration/HumanColonCancerPatient1/cell_metadata.csv'
metadata = read.csv(path2, row.names = 1)

location = metadata[,c('center_x','center_y')]
location2 = as.matrix(location)

good_genes = (colnames(count2) %in% colnames(count1))

count2 = count2[,good_genes]

idx = rowSums(count2) >= 20
write.csv(idx,file='data/cellidx_patient1.csv',row.names=F)

location2 = location2[order(as.numeric(rownames(metadata))),]
location2 = location2[idx,]

count2 = count2[idx,]

batch_id2 = groupondist(location2,size = 5000, lengthscale = 100)
batch_id2 = paste0("2_",batch_id2)

# load dataset 3
path1 = '/gladstone/engelhardt/pelka-collaboration/HumanColonCancerPatient2/cell_by_gene.csv'
count = read.csv(path1, row.names = 1)
count3 = as.matrix(count)
path2 = '/gladstone/engelhardt/pelka-collaboration/HumanColonCancerPatient2/cell_metadata.csv'
metadata = read.csv(path2, row.names = 1)

location = metadata[,c('center_x','center_y')]
location3 = as.matrix(location)

good_genes = (colnames(count3) %in% colnames(count1))

count3 = count3[,good_genes]

idx = rowSums(count3) >= 20
write.csv(idx,file='data/cellidx_patient2.csv',row.names=F)

location3 = location3[order(as.numeric(rownames(metadata))),]
location3 = location3[idx,]
count3 = count3[idx,]

batch_id3 = groupondist(location3,size = 5000, lengthscale = 100)
batch_id3 = paste0("3_",batch_id3)

sum(colnames(count1) != colnames(count2))
sum(colnames(count1) != colnames(count3))

count2 = count2[,colnames(count1)]
count3 = count3[,colnames(count1)]

# all together
total_batch = c(batch_id1,batch_id2,batch_id3)
total_count = rbind(count1,count2,count3)
total_location = rbind(location1,location2*9.259,location3*9.259)

# start = Sys.time()
# out = nmfspatial_batch(total_count, 5, location = total_location, batch = total_batch, lengthscale = 1000, initial = 1, maxiter = 1000, tolerance = 1e-6, error_freq = 10, normalize = FALSE)
# stop = Sys.time()
# out$time = stop - start
# print(stop - start)
# save(out, file = paste0("modelssaved/three_sample_f20_s5K_l500_i1000.RData"))

for(s in c(50,200,500,1000,5000)){
print(s)
start = Sys.time()
out = nmfspatial_batch(total_count, 20, location = total_location, batch = total_batch, lengthscale = s, initial = 1, maxiter = 1000, tolerance = 1e-6, error_freq = 10, kernel_cutoff = 0.1)
stop = Sys.time()
out$time = stop - start
print(stop - start)

save(out, file = paste0("modelssaved/three_sample_f20_s5K_l",s,"_norm_co02.RData"))
}






# row_sum = rowSums(total_count)
# total_count = total_count/row_sum*mean(row_sum) 

# start = Sys.time()
# out = nmfgen(total_count, 20, iter = 1000)
# stop = Sys.time()
# out$time = stop - start
# print(stop - start)
# save(out, file = paste0("modelssaved/three_sample_f20_nmfgen_norm_iter1000.RData"))


