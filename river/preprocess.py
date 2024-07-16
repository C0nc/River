# uncomment this if you want to use interactive plot (only works in Jupyter not works in VScode)
# %matplotlib widget

import scanpy as sc
import numpy as np
import pandas as pd

import scSLAT
from scSLAT.model import Cal_Spatial_Net, load_anndatas, run_SLAT, spatial_match
from scSLAT.viz import match_3D_multi, hist, Sankey
from scSLAT.metrics import region_statistics\
    
#adata1 = adata_list[33]#adata_list[1]
#adata2 = adata_list[0]

def seed_everything(seed: int):
    import random, os
    import numpy as np
    import torch
    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


def do_slat_pair(adata1, adata2, feature='pca'):

    Cal_Spatial_Net(adata1, k_cutoff=20, model='KNN')
    Cal_Spatial_Net(adata2, k_cutoff=20, model='KNN')

    edges, features = load_anndatas([adata1, adata2], feature=feature)

    embd0, embd1, time = run_SLAT(features, edges)

    best, index, distance = spatial_match(features, adatas=[adata1, adata2], reorder=False)
    
    return [adata1, adata2], [best]




def get_the_feature(adata_list, matching_list):
    
    overlap = set(adata_list[0].var_names.tolist())
    
    for adata in adata_list[1:]:
        
        overlap &= set(adata.var_names.tolist())
        
    overlap = sorted(list(overlap))
    
    print(len(overlap))

    X1 = adata_list[0][:,overlap].X.toarray()

    spatial = adata_list[0][:,overlap].obsm['spatial']
        
    X = [X1]
    
    spatials = [spatial]
    
    y = [np.array([0] * adata_list[0].X.shape[0])]
    
    
    
    for i, adata in enumerate(adata_list[1:]):
        
        if len(matching_list[i].shape)>1:
        
            spatials.append(spatial[matching_list[i][1]])
            
            X.append(adata[matching_list[i][0],overlap].X.toarray())
            
            y.append(np.array([i+1] * adata[matching_list[i][0],overlap].X.shape[0]))
        
        else:

            
            spatials.append(spatial[matching_list[i]])
            
            X.append(adata[:,overlap].X.toarray())
            
            y.append(np.array([i+1] * adata[:,overlap].X.shape[0]))

    
    return np.concatenate(X, axis=0), np.concatenate(spatials, axis=0), np.concatenate(y, axis=0).astype(int), overlap 


def get_the_multi_feature(adata_list, matching_list):
    
    overlap = set(adata_list[0].var_names.tolist())
    
    for adata in adata_list[1:]:
        
        overlap &= set(adata.var_names.tolist())
        
    overlap = sorted(list(overlap))

    X1 = adata_list[0][:,overlap].X.toarray()

    y1 = adata_list[0][:,overlap].obsm['spatial']
    
    X1_addition = np.ones((y1.shape[0], 1)) * 0
    
    X = [np.concatenate((X1, X1_addition), axis=1)]
    
    y = [y1]
    
    y_1 = [np.array([0] * adata_list[0].X.shape[0])]
    
    
    for i, adata in enumerate(adata_list[1:]):
        
        if len(matching_list[i].shape)>1:
        
            tmp_y1 = y1[matching_list[i][1]]
        
            tmp_y1[...,-1] = np.ones((tmp_y1.shape[0])) * (i+1)
        
            y.append(tmp_y1)
            
            X.append(adata[matching_list[i][0],overlap].X.toarray())
            
            y_1.append(np.array([i+1] * adata[matching_list[i][0],overlap].X.shape[0]))
        
        else:
            
            tmp_y1 = y1[matching_list[i]]
            
            X1_addition = np.ones((tmp_y1.shape[0], 1)) * (i+1)
            
            y.append(tmp_y1)
            
            X.append(np.concatenate((adata[:,overlap].X.toarray(), X1_addition), axis=1))
            
            y_1.append(np.array([i+1] * adata[:,overlap].X.shape[0]))
  
    return  np.concatenate(X, axis=0), np.concatenate(spatials, axis=0), np.concatenate(y, axis=0).astype(int), adata_list[0].var_names.tolist() 


def get_the_sim_feature(adata_list):
    
    
    X = [adata_list[0].X]
    
    y1 = adata_list[0].obsm['spatial']

    y = [y1]

    y_1 = [np.array([0] * adata_list[0].X.shape[0])]

    for i, adata in enumerate(adata_list[1:]):
        
        y.append(adata.obsm['spatial'])
        
        X.append(adata.X)
        
        y_1.append(np.array([i+1] * adata.X.shape[0]))

    return np.concatenate(X, axis=0), np.concatenate(y, axis=0), np.concatenate(y_1, axis=0).astype(int), adata_list[0].var_names.tolist() 
