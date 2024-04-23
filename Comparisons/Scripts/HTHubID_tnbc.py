#%%
# import packages and set paths
import argparse
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
#import geopandas as gpd
import random
import seaborn as sns
import os
import sys
os.chdir("~/")
sample_name = "tnbc_t_cells_immune_neighborhood_noshuffle"
assay = 'MIBITOF'

## preprocess data
import pickle
df_dict=pickle.load(open("./github/spatial_lda/data/tnbc//patient_dfs.pkl", "rb"))


ss=0
for ii in list(df_dict.keys()):
    df_dict[ii].index=[(ii, df_dict[ii].index[i]) for i in range(0, len( df_dict[ii]))]

# Extract the values (DataFrames) from the dictionary
dfs = list(df_dict.values())

# Concatenate the DataFrames along the rows (axis=0)
df0 = pd.concat(dfs)

df=df0.iloc[:,2:44]
# print(df_dict[1].columns)
# print(cells_features.columns)

print(df)


#mtx=cell_feature.copy()
mtx=df.copy()
#dat=np.array(cell_feature)
coordinates=np.array(df0[['x','y']])
coordinates

#%%
# flip coordinate
cor1=coordinates[:,1].copy()
cor0=coordinates[:,0].copy()
coordinates[:,1]=cor0
coordinates[:,0]=cor1
coordinates


# %%
# check average count
sum(mtx.sum(axis=1)<1)

# %%
mtx
index_new=[]
for ii in mtx.index:
    index_new.append(str(ii[0])+"_"+str(ii[1]))

mtx.index=index_new
#coordinates.index=index_new
print(index_new)


## shuffle the rows
mtx0=mtx.copy()
coordinates0=coordinates.copy()

np.random.seed(2024)
#shuffled_index = mtx.index.to_numpy()
shuffled_index=np.random.permutation(range(mtx.shape[0]))
#np.random.shuffle(shuffled_index)
mtx = mtx.iloc[shuffled_index,:]
coordinates = coordinates[shuffled_index]

# %%
mtx.columns

adata = sc.AnnData(mtx, 
    mtx.index.to_frame(), 
    mtx.columns.to_frame())

adata.obsm['spatial']=coordinates
df0.index=index_new
df0=df0.iloc[shuffled_index,:]

adata.obs=pd.concat([adata.obs, df0], axis=1)
adata.obs['istumor']=~adata.obs['isimmune']
gene_names=adata.var.iloc[:,0].tolist()
print(gene_names[0:5])
#%5

# %%
## HTHubID
seed = 2024
min_genes = 0
min_counts = 3

outpath = './'
# %%
adata_sub=adata#sc.pp.subsample(adata, fraction=0.1, copy=True, random_state=0)
if min_genes>0:
    sc.pp.filter_cells(adata_sub, min_genes=min_genes)
elif min_counts>0:
    sc.pp.filter_cells(adata_sub,  min_counts= min_counts)

#%%
# save the object
if outpath == './':
    os.makedirs(outpath+"/"+str(sample_name)+"/Tables/", exist_ok=True)
    os.makedirs(outpath+"/"+str(sample_name)+"/Plots/", exist_ok=True)
    os.makedirs(outpath+"/"+str(sample_name)+"/Objects/", exist_ok=True)
    outpath = outpath+"/"+str(sample_name)+"/Objects/" + str(sample_name)+'_'+str(assay)+"_adata_mingenes"+str(min_genes)+"_mincounts"+str(min_counts)+".h5ad"

adata_sub.obs.rename(columns={0: "index"}, inplace=True)
adata_sub.var.rename(columns={0: "gene"}, inplace=True)
adata_sub.write_h5ad(outpath)


#%%


#### RUN HTHubID
from HTHubID.EstCellPrograms import *
import os   
# min_genes = 0
# min_counts = 0


#anndata_path="/gladstone/engelhardt/home/chwu/Research/HubID/Gensim/Objects/"+str(sample_name)+"_adata_"+str(assay)+".h5ad"
anndata_path='./'+str(sample_name)+'/Objects/'+str(sample_name)+'_'+str(assay)+'_adata_mingenes'+str(min_genes)+"_mincounts"+str(min_counts)+".h5ad"
outdir = './'
random_state = 2024
fraction = 1
ncol = 5
spot_size_cluster = 20
spot_size_activity = 20
is_filter=True


# import files
adata_sub = sc.read_h5ad(anndata_path)

# filter only tumor cells
#adata_sub=adata_sub[adata_sub.obs['isimmune']==False]


#%%
# run EstLandscape
import time
t0 = time.time()
cellprogram_estimator = CellProgramEstimator(adata_sub=adata_sub, sample_name=sample_name, assay=assay, outdir=outdir, is_filter=is_filter, random_state=random_state, fraction=fraction, ncol=ncol, spot_size_cluster=spot_size_cluster, spot_size_activity=spot_size_activity, multi_samples=True, plot_set='4')
# Call the EstLandscape method
cellprogram_estimator.EstCellPrograms()
t1 = time.time()
total=t1-t0


##################################
###### Landscape estimation#########
# landscape estimation
# importpackage
import os
from HTHubID.EstLandscapes import *
os.chdir("/gladstone/engelhardt/lab/cywu/HubID/SpatialLDA")

# sample_name = "tnbc"
# assay = 'MIBITOF'

# from last step
min_genes = 0
min_counts = 3

## double check function for baysor
# importpackage
anndata_path='./'+str(sample_name)+'/Objects/'+str(sample_name)+'_'+str(assay)+'_adata_mingenes'+str(min_genes)+"_mincounts"+str(min_counts)+".h5ad"

#prop_path='./'+str(sample_name)+'/Tables/prop_'+str(sample_name)+'_'+str(assay)+'.csv'
#program_path='./'+str(sample_name)+'/Tables/topic_info_weighted_'+str(sample_name)+'_'+str(assay)+'.csv'
neighbor_mode='radius'
#neighbor_mode='NNeighbors'
eps=100
col_type='program'
outdir='./'
n_components=5
#assay='merFISH_cellpose'
init='random'
random_state=2024#0
alpha_W=0 
fraction=1
ncol=5
spot_size= 20 #70#600 #50  
n_neighbors=0

# import files
adata_sub = sc.read_h5ad(anndata_path)
prop=adata_sub.X
cellid=adata_sub.obs.index
#prop=prop.iloc[:,1:]
topicid=df.columns#prop.columns
#prop=np.array(prop)
#prop=adata_sub.X


# remove NA
rm_id=np.union1d(np.where((np.isnan(prop).any(axis=1))), np.where(pd.DataFrame(adata_sub.obsm['spatial']).isna().any(axis=1)))
sel_id=np.delete(np.array(range(prop.shape[0])),rm_id)
adata_sub=adata_sub[sel_id,:]
prop=prop[sel_id,:]
cellid=cellid[sel_id]

#prop=adata_sub.X


## set the colors
#assay='MIBITOF_raw_values_for_est'
#color_map= {'NN': 'grey'}
#custom_cmap = plt.cm.colors.ListedColormap([color_map.get(category, 'default') for category in adata.obs['category'].unique()])

# Define colors for a specific group manually
manual_color = 'grey'

# Get unique groups from your data
groups= [str(i) for i in range(5)]
unique_groups = np.array(groups)

# Generate the Paired color map
paired_cmap = plt.cm.get_cmap('Set1')

# Create a dictionary to store the palette
full_palette = {}

# Assign manual color to the specified group
full_palette['NN'] = manual_color

# Assign Paired colormap to other groups
for i, group in enumerate(unique_groups[0:], start=1):
    full_palette[group] = paired_cmap(i)

# Plot spatial data with the full palette
#sc.pl.spatial(adata, color=group_column, palette=full_palette)



# run EstLandscape
#landscape_estimator = LandscapeEstimator(adata_sub=adata_sub, prop=prop, program=program, sample_name=sample_name, assay=assay, cellid=cellid, neighbor_mode=neighbor_mode, n_neighbors=n_neighbors, eps=eps, outdir=outdir, n_components=n_components, init=init, random_state=random_state,alpha_W=alpha_W, fraction=fraction, ncol=ncol, spot_size=spot_size)
#landscape_estimator = LandscapeEstimator(adata_sub=adata_sub, prop=prop, program=program, sample_name=sample_name, assay=assay, cellid=cellid, topicid=topicid, neighbor_mode=neighbor_mode, n_neighbors=n_neighbors, eps=eps, col_type=col_type, outdir=outdir, n_components=n_components, init=init, random_state=random_state,alpha_W=alpha_W, fraction=fraction, ncol=ncol, spot_size=spot_size, multi_samples=True, tumor_only_est=True)# tumor only just for landscape est not for smoothing
landscape_estimator = LandscapeEstimator(adata_sub=adata_sub, prop=prop, sample_name=sample_name, assay=assay, cellid=cellid, topicid=topicid, neighbor_mode=neighbor_mode, n_neighbors=n_neighbors, eps=eps, col_type=col_type, outdir=outdir, n_components=n_components, init=init, random_state=random_state,alpha_W=alpha_W, fraction=fraction, ncol=ncol, spot_size=spot_size, multi_samples=True, tumor_only_est=True, palette=full_palette)#
# Call the EstLandscape method
landscape_estimator.EstLandscape()





