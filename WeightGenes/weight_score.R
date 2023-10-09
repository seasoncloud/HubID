#score matrix
dat=read.csv("./topic_info_hdp_LH_VMSC02001_all_seed0.csv", stringsAsFactors = F, row.names = 1) 
dat=t(dat) # we want gene by topic


## remove genes that are all zero
gene_prop=rowSums(dat)
sel_idx=which(gene_prop!=0)

dat_new=NULL
for(ii in 1:nrow(dat[sel_idx,,drop=F])){
  rr=dat[sel_idx,,drop=F][ii,]
  m1=max(rr)
  m2=max(rr[-which(rr==m1)])
  mm=rep(m1, length(rr))
  mm[which(rr==m1)]=m2
  ns=rr*log((rr+1)/(mm+1))
  dat_new=rbind(dat_new, ns)
  print(ii)
}
rownames(dat_new)=rownames(dat)[sel_idx]

dat_new=apply(dat_new, c(1,2), function(x) round(x,10))

write.csv(dat_new, "./topic_info_all_seed0_weighted_rm0.csv", quote = F, row.names =T )


