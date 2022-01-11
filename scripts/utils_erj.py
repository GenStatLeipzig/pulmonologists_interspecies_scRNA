import scanpy as sc
import scvelo as scv
import pandas as pd
import numpy as np


def calculate_velocity(path_to_loom):
    '''
    parameters:
        path_to_loom: string, path to loom file (e.g.'loom_lung_1.loom')
        
    output:
        adata: anndata object with spliced and unspliced entries
    '''
    adata = scv.read(path_to_loom,cache = True)
    adata.var_names_make_unique()
    scv.pp.moments(adata)
    scv.pp.neighbors(adata)
    scv.tl.velocity(adata, mode='stochastic')
    scv.tl.velocity_graph(adata)
    scv.pl.velocity_embedding_stream(adata)
    return adata

def filter_adata_obs(adata,obs_feature,obs_value):
    '''
    parameters:
        adata         : anndata object
        obs_feature   : string, column name in adata.obs
        obs_value     : string or double or float or int, value in adata.obs['obs_feature'] to filter for
   
   output:
       adata_filtered : anndata object
    '''
    adata_filtered = adata[adata.obs[obs_feature] == obs_value]
    return adata_filtered

def create_ID_velo(adata,start_pos_base):
    '''
    parameters:
        adata            : anndata object
        start_pos_base   : int, position of base sequence in adata.obs.index  
    '''
    adata.obs['ID_velo'] = adata.obs.index
    for j in range(len(adata)):
        adata.obs['ID_velo'][j] = adata.obs['ID_velo'][j][start_pos_base:start_pos_base+16]

def find_cell_tags(adata_velo,adata_erj):
    '''
    parameters:
        adata_velo  : anndata object, output of loom file
        adata_erj   : anndata object, already preprocessed 
        
    return:
        used_index  : list of indices of adata_velo that are also used in preprocessed version
    '''
    used_index = []
    used_tag = list(adata_erj.obs['ID_velo'].values)
    for q in range(len(adata_velo)):
        if adata_velo.obs['ID_velo'][q] in used_tag:
            used_index.append(q)
    return used_index

def merge_obsm(adata_velo,adata_erj,obsm_string):
    '''
    parameters:
        adata_velo  : anndata object, output of loom file
        adata_erj   : anndata object, already preprocessed 
        obsm_string : string, column name of adata_erj.obsm
    '''
    pos_dict_velo = {}
    for p in range(len(adata_velo)):
         pos_dict_velo[adata_velo.obs['ID_velo'][p]] = p
    obsm_object = adata_erj.obsm[obsm_string].copy()
    for q in range(len(adata_erj)):
         pos_new = pos_dict_velo[adata_erj.obs['ID_velo'].values[q]]
         obsm_object[pos_new,:] = adata_erj.obsm[obsm_string][q,:].copy()
    adata_velo.obsm[obsm_string] = obsm_object.copy()

def velocity_analysis_erj(path_to_loom,
                         path_to_h5ad,
                         file_ident,
                         char_to_skip_loom = 12,
                         char_to_skip_h5ad = 13,
                         col_name_obs = 'orig.ident',
                         col_value_obs = 'ma_d0_lung_1',
                         write_h5ad = True,
                         save_velo_graphic = True):
    '''
    parameters:
        path_to_loom      : string
        path_to_h5ad      : string
        file_ident        : string, prefix(name to save)
        char_to_skip_loom : int, nr of characters in obs.index (loom) until base sequence (e.g AAATGC..)
        char_to_skip_h5ad : int, nr of characters in obs.index (h5ad) until base sequence (e.g AAATGC..)
        write_h5ad        : True or False, whether to save final anndata object as .h5ad
        save_velo_graphic : True or False, whether to save RNA velocity graphics
    '''
    scv.set_figure_params()
    h1 = calculate_velocity(path_to_loom)
    erj = sc.read_h5ad(path_to_h5ad)
    erj_2 = filter_adata_obs(erj,col_name_obs,col_value_obs)
    create_ID_velo(erj_2,char_to_skip_h5ad)
    create_ID_velo(h1,char_to_skip_loom)
    ind = find_cell_tags(h1,erj_2)
    adata_fin = h1[ind]
    adata_fin.obs = adata_fin.obs.merge(erj_2.obs[['species_tissue','tissue_v1','cluster_seurat_v1',
                                           'populations','spezies','ID_velo']], on = 'ID_velo')
    merge_obsm(adata_fin,erj_2,'X_umap')
    if write_h5ad == True:
        adata_fin.write(file_ident+'_velo.h5ad')
    if save_velo_graphic == True:
        scv.pl.velocity_embedding_stream(adata_fin,
                                         basis='X_umap',
                                         color = 'populations',
                                         linewidth = 1.4,
                                         save=file_ident+'_velo_populations.png')
        scv.pl.velocity_embedding_stream(adata_fin, basis='X_umap',
                                         color = 'cluster_seurat_v1',
                                         linewidth = 1.4,
                                        save=file_ident+'_velo_seurat_cluster.png')
    else:
        scv.pl.velocity_embedding_stream(adata_fin,
                                         basis='X_umap',
                                         color = 'populations',
                                         linewidth = 1.4)
        scv.pl.velocity_embedding_stream(adata_fin, basis='X_umap',
                                         color = 'cluster_seurat_v1',
                                         linewidth = 1.4)