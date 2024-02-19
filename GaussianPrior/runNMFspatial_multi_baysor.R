
source("model/NMFbatch.R")
batch_size = 10000

# load dataset 1
path1 = '/gladstone/engelhardt/pelka-collaboration/HuColonCa-FFPE-ImmuOnco-LH_VMSC02001_20220427/baysor_20231207/baysor_res_merged/cell_by_gene.csv.gz'
count = read.csv(gzfile(path1))
count1 = as.matrix(count)
genes = colnames(count)
path1 = '/gladstone/engelhardt/pelka-collaboration/HuColonCa-FFPE-ImmuOnco-LH_VMSC02001_20220427/baysor_20231207/baysor_res_merged/segmentation_cell_stats.csv.gz'
metadata = read.csv(gzfile(path1))

location = cbind(metadata$x, metadata$y)
location1 = as.matrix(location)

idx = !(rowSums(count1) < 20 | is.na(location1[,1]))

write.csv(idx,file='data/cellidx_baysor1.csv',row.names=F)

location1 = location1[idx,]

count1 = count1[idx,]

batch_id1 = groupondist(location1,size = batch_size)
batch_id1 = paste0("1_",batch_id1)

# load dataset 2
path1 = '/gladstone/engelhardt/pelka-collaboration/HumanColonCancerPatient1/cell_by_gene.csv'
count = read.csv(path1, row.names = 1)
count2 = as.matrix(count)
path2 = '/gladstone/engelhardt/pelka-collaboration/HumanColonCancerPatient1/cell_metadata.csv'
metadata = read.csv(path2, row.names = 1)

location = metadata[,c('center_x','center_y')]
location2 = as.matrix(location)

count2 = count2[,genes]

idx = rowSums(count2) >= 20
write.csv(idx,file='data/cellidx_patient1.csv',row.names=F)

location2 = location2[order(as.numeric(rownames(metadata))),]
location2 = location2[idx,]

count2 = count2[idx,]

batch_id2 = groupondist(location2,size = batch_size)
batch_id2 = paste0("2_",batch_id2)

# load dataset 3
path1 = '/gladstone/engelhardt/pelka-collaboration/HumanColonCancerPatient2/cell_by_gene.csv'
count = read.csv(path1, row.names = 1)
count3 = as.matrix(count)
path2 = '/gladstone/engelhardt/pelka-collaboration/HumanColonCancerPatient2/cell_metadata.csv'
metadata = read.csv(path2, row.names = 1)

location = metadata[,c('center_x','center_y')]
location3 = as.matrix(location)

count3 = count3[,genes]

idx = rowSums(count3) >= 20
write.csv(idx,file='data/cellidx_patient2.csv',row.names=F)

location3 = location3[order(as.numeric(rownames(metadata))),]
location3 = location3[idx,]
count3 = count3[idx,]

batch_id3 = groupondist(location3,size = batch_size)
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
# out = nmfspatial_batch(total_count, 20, location = total_location, batch = total_batch, lengthscale = 5000, initial = 1, maxiter = 1000, tolerance = 1e-6, error_freq = 10, normalize = TRUE, kernel_cutoff = 0.1)
# stop = Sys.time()
# out$time = stop - start
# print(stop - start)
# save(out, file = paste0("modelssaved/three_sample_f20_s10K_l500_baysor_co01.RData"))
start = Sys.time()
out = nmfspatial_batch(total_count, 20, location = total_location, batch = total_batch, lengthscale = 100, initial = 1, maxiter = 1000, tolerance = 1e-6, error_freq = 10, kernel_cutoff = 0.5)
stop = Sys.time()
out$time = stop - start
print(stop - start)

save(out, file = paste0("modelssaved/three_sample_f20_s10K_l100_norm_baysor_co05_v2.RData"))


for(s in c(0.1,0.5)){
print(s)
start = Sys.time()
out = nmfspatial_batch(total_count, 20, location = total_location, batch = total_batch, lengthscale = 50, initial = 1, maxiter = 1000, tolerance = 1e-6, error_freq = 10, kernel_cutoff = s)
stop = Sys.time()
out$time = stop - start
print(stop - start)

save(out, file = paste0("modelssaved/three_sample_f20_s10K_l50_norm_baysor_co",s,".RData"))
}


# row_sum = rowSums(total_count)
# total_count = total_count/row_sum*mean(row_sum) 

# start = Sys.time()
# out = nmfgen(total_count, 20, iter = 1000)
# stop = Sys.time()
# out$time = stop - start
# print(stop - start)
# save(out, file = paste0("modelssaved/three_sample_f20_nmfgen_norm_iter1000.RData"))


