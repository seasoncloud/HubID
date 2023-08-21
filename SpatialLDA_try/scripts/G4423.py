import os
import sys

os.chdir("/wynton/home/engelhardt/chwu/Projects/HubID/SpatialLDA")
PATH_TO_MODELS = f'./models/'
#PATH_TO_DF_PKL = f'./data//spleen_df.pkl'
PATH_TO_FEATURES_PKL = f'./data/cells_features.pkl'

# import packages
%load_ext autoreload
%autoreload 2

import functools
import logging
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import pickle
import scipy
import seaborn as sns
from sklearn.model_selection import train_test_split
import time
import tqdm
# Spatial LDA imports
from spatial_lda.featurization import neighborhood_to_cluster, neighborhood_to_marker
from spatial_lda.featurization import make_nearest_neighbor_graph
from spatial_lda.featurization import make_merged_difference_matrices
from spatial_lda.featurization import featurize_samples
from spatial_lda.visualization import plot_samples_in_a_row
from spatial_lda.visualization import plot_bcell_topic_multicolor
import spatial_lda.model

N_PARALLEL_PROCESSES = 8#@param
TRAIN_SIZE_FRACTION = 0.99 #@param
N_TOPICS_LIST = [3, 5, 8, 10] #@param

logger = logging.getLogger()
logger.setLevel(logging.INFO)

# visualization parameters
blue = sns.color_palette()[0]
green = sns.color_palette()[1]
red = sns.color_palette()[2]
sns.set_context("paper",font_scale=1)
sns.set_style('white')
sns.set_context("notebook", font_scale=1.0, rc={"lines.linewidth": 1.0})


# import data
mtx=scipy.sparse.load_npz("../neighbor_smooth/Data/G4423_15_neighbor_sums.npz")
with open("../neighbor_smooth//Data/G4423_15_neighbor_sums_colnames.csv", "r+") as file1:
    # Reading form a file
    colnames=file1.read()
    
colnames=colnames.strip("\n").split(",")

mtx=pd.DataFrame(mtx.todense())
mtx.columns=colnames
mtx=mtx.iloc[:,1:]

# import spatial info
meta=pd.read_csv('/wynton/group/gladstone/users/cywu/Broad_MERFISH_Datasets/G4423_20220427/cellpose_cell_metadata.csv')
coordinates=np.array(meta[['center_x', 'center_y']])
# flip y
coordinates[:,1]=max(coordinates[:,1])-coordinates[:,1]

mtx['isb']=True
mtx['sample']='G4423'
mtx['sample.X']=coordinates[:,0]
mtx['sample.Y']=coordinates[:,1]

# subsample and generate dic
import random
idx=random.sample(range((mtx.shape[0])), 10000)

df_dict = {'G4423': mtx.loc[idx]}


# Featurize the data
#%%time
for df in df_dict.values():
    df['x'] = df['sample.X']
    df['y'] = df['sample.Y']
#wt_samples = [ x for x in codex_df_dict.keys() if x.startswith("BALBc")]
#spleen_dfs = dict(zip(wt_samples, [ codex_df_dict[x] for x in wt_samples]))
  
cells_features = featurize_samples(df_dict, neighborhood_to_marker, 150, 'isb',
                             'sample.X', 'sample.Y', include_anchors=True,
                             n_processes=N_PARALLEL_PROCESSES)

#featurize_spleens(spleen_dfs, neighborhood_to_cluster, radius=100,n_processes=N_PARALLEL_PROCESSES)
with open(PATH_TO_FEATURES_PKL, 'wb') as f:
    pickle.dump(cells_features, f)

# compute difference matrices
difference_matrices = make_merged_difference_matrices(cells_features, df_dict,
                                                             'sample.X', 'sample.Y')
all_sample_idxs = cells_features.index.map(lambda x: x[0])
_sets = train_test_split(cells_features, 
                         test_size=1. - TRAIN_SIZE_FRACTION,
                         stratify=all_sample_idxs)
train_cells_features, test_cells__features = _sets
train_difference_matrices = make_merged_difference_matrices(
    train_cells_features, df_dict,
    'sample.X', 'sample.Y')
train_idxs = train_cells_features.index.map(lambda x: x[0])


## visualization
from spatial_lda.visualization import plot_adjacency_graph

def make_plot_fn(difference_matrices):  
    def plot_fn(ax, tumor_idx, features_df, patient_dfs):
        plot_adjacency_graph(ax, tumor_idx, features_df, patient_dfs, difference_matrices)
    return plot_fn
_plot_fn = make_plot_fn(difference_matrices)

# plot_samples_in_a_row(cells_features, _plot_fn, df_dict)

from spatial_lda.visualization import plot_samples_in_a_row
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
from matplotlib.backends.backend_pdf import PdfPages
from spatial_lda.visualization import plot_samples_in_a_row
with PdfPages('./plots/G4423_all_adjacency.pdf') as pdf_pages:
    # df1 = dftest.select_dtypes([np.int, np.float, np.object])
    # for i, col in enumerate(N_TOPICS_LIST):
    #     figu = plt.figure(i)
    #     plot = sns.countplot(x=col, data=df1)
    #     pdf_pages.savefig(figu)
    plt.figure()
    plot_samples_in_a_row(cells_features, _plot_fn, df_dict)
    pdf_pages.savefig()
    plt.close()

## Spatial LDA results
N_TOPICS_LIST

from spatial_lda.model import order_topics_consistently
spatial_lda_models = {}  
difference_penalty = 0.25  
for n_topics in N_TOPICS_LIST:
  path_to_train_model = '_'.join((f'{PATH_TO_MODELS}/training_all',
                                  f'penalty={difference_penalty}',
                                  f'topics={n_topics}',
                                  f'trainfrac={TRAIN_SIZE_FRACTION}')) + '.pkl'
  if not os.path.exists(path_to_train_model):
    print(f'Running n_topics={n_topics}, d={difference_penalty}\n')
    spatial_lda_model = spatial_lda.model.train(sample_features=train_cells_features, 
                                                difference_matrices=train_difference_matrices,
                                                difference_penalty=difference_penalty,
                                                n_topics=n_topics,
                                                n_parallel_processes=N_PARALLEL_PROCESSES,                                                                         
                                                verbosity=1,
                                                admm_rho=0.1,
                                                primal_dual_mu=1e+5)
    spatial_lda_models[n_topics] = spatial_lda_model
    with open(path_to_train_model, 'wb') as f:
      pickle.dump(spatial_lda_model, f)    
  else:
    with open(path_to_train_model, 'rb') as f:
      spatial_lda_models[n_topics] = pickle.load(f)
      
order_topics_consistently(spatial_lda_models.values())     

lda_3 = spatial_lda_models[3]
topic_weights_3 = lda_3.topic_weights
lda_5 = spatial_lda_models[5]
topic_weights_5 = lda_5.topic_weights
lda_8 = spatial_lda_models[8]
topic_weights_8 = lda_8.topic_weights
lda_10 = spatial_lda_models[10]
topic_weights_10 = lda_10.topic_weights
samples = ['G4423']


## visualization
from spatial_lda.visualization import plot_samples_in_a_row
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
from matplotlib.backends.backend_pdf import PdfPages
from spatial_lda.visualization import plot_samples_in_a_row
with PdfPages('./plots/G4423_all_topic_weights.pdf') as pdf_pages:
    # df1 = dftest.select_dtypes([np.int, np.float, np.object])
    # for i, col in enumerate(N_TOPICS_LIST):
    #     figu = plt.figure(i)
    #     plot = sns.countplot(x=col, data=df1)
    #     pdf_pages.savefig(figu)
    plt.figure()
    plot_samples_in_a_row(topic_weights_3, plot_bcell_topic_multicolor, df_dict, tumor_set=samples)
    pdf_pages.savefig()
    plt.close()
    
    plt.figure()
    plot_samples_in_a_row(topic_weights_5, plot_bcell_topic_multicolor, df_dict, tumor_set=samples)
    pdf_pages.savefig()
    plt.close()
    
    plt.figure()
    plot_samples_in_a_row(topic_weights_8, plot_bcell_topic_multicolor, df_dict, tumor_set=samples)
    pdf_pages.savefig()
    plt.close()
    
    plt.figure()
    plot_samples_in_a_row(topic_weights_10, plot_bcell_topic_multicolor, df_dict, tumor_set=samples)
    pdf_pages.savefig()
    plt.close()
    
    



