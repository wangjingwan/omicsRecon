# conda activate cell2loc_env
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import cell2location
import scvi
import sys
from scipy.sparse import csr_matrix
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-s', '--st-file', dest='st_file', required=True, help='Spatial transcriptomics data')
parser.add_argument('-c', '--st-coordinate', dest='st_coord', required=False, help='Spatial coordinates of the spatial transcriptomics data')
parser.add_argument('-r', '--sc-file', dest='sc_file', required=True, help='Single-cell candidate library of the corresponding ST tissue')
parser.add_argument('-m', '--meta-file', dest='meta_file', required=True, help='Cell-type annotation of the single-cell candidate library')
parser.add_argument('-t', '--type-key', dest='tp_key', required=True, help='The colname of celltype in SC_meta')
parser.add_argument('-o', '--out-dir', dest='out_dir', required=False, help='Output file path')
parser.add_argument('-b', '--batch-key', dest='batch_key', required=False, help='batch column name in SC meta')

args = parser.parse_args()
def read_csv_tsv(filename):
    if 'csv' in filename:
        tmp = pd.read_csv(filename, sep = ',',header = 0,index_col=0)
    else:
        tmp = pd.read_csv(filename, sep = '\t',header = 0,index_col=0)
    return tmp

# 1. results_folder
if args.out_dir is not None:
    results_folder = args.out_dir + '/' 
    if not os.path.exists(results_folder):
        os.makedirs(results_folder)
else:
    results_folder = os.getcwd() + '/'
print(f'Saving the c2l results in {results_folder}')
ref_run_name = f'{results_folder}/reference_signatures'
run_name = f'{results_folder}/cell2location_map'
# 1.ST
st_mat = read_csv_tsv(args.st_file)
st_meta = read_csv_tsv(args.st_coord)
if 'row' in st_meta.columns.tolist():
    # if the x and y coordinates are indexed with row and col
    st_meta = st_meta[['row','col']].copy()
else:
    st_meta = st_meta[['x','y']].copy()
st_meta.columns = ['x','y']
st_exp = anndata.AnnData(csr_matrix(st_mat))
st_exp.obs = st_meta
st_exp.var = pd.DataFrame(st_mat.columns, columns=['symbol'])
st_exp.var_names = st_exp.var['symbol'].copy()
sc.pp.filter_cells(st_exp, min_genes=200)
sc.pp.filter_genes(st_exp, min_cells=3)
st_exp.obs['sample'] = 0
# 2.SC
sc_mat = read_csv_tsv(args.sc_file)
meta_df = read_csv_tsv(args.meta_file)
sc_mat = sc_mat.loc[meta_df.index]
# meta_df.index = sc_mat.index
sc_mat = np.round(sc_mat,0)
adata_sc = anndata.AnnData(csr_matrix(sc_mat))
adata_sc.obs = meta_df
adata_sc.var = pd.DataFrame(sc_mat.columns, columns=['symbol'])
adata_sc.var_names = adata_sc.var['symbol'].copy()
sc.pp.filter_cells(adata_sc, min_genes=200)
sc.pp.filter_genes(adata_sc, min_cells=3)
# 
if args.batch_key is not None:
    batch_key = args.batch_key
else:
    adata_sc.obs['sample'] = 0
    batch_key = 'sample'
######sc_loc = f'{inputDir}/SC/sc_adata.h5ad'
from cell2location.utils.filtering import filter_genes
selected = filter_genes(adata_sc, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
# filter the object
adata_sc = adata_sc[:, selected].copy()
scvi.data.setup_anndata(adata=adata_sc,
                        # 10X reaction / sample / batch
                        batch_key=batch_key,
                        # cell type, covariate used for constructing signatures
                        labels_key=args.tp_key
                       )
scvi.data.view_anndata_setup(adata_sc)

from cell2location.models import RegressionModel
mod = RegressionModel(adata_sc)
# Use all data for training (validation not implemented yet, train_size=1)
mod.train(max_epochs=250, batch_size=2500, train_size=1, lr=0.002, use_gpu=True)
# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_sc = mod.export_posterior(
    adata_sc, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}
)
mod.save(f"{ref_run_name}", overwrite=True)
if 'means_per_cluster_mu_fg' in adata_sc.varm.keys():
    inf_aver = adata_sc.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}' for i in adata_sc.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_sc.var[[f'means_per_cluster_mu_fg_{i}' for i in adata_sc.uns['mod']['factor_names']]].copy()

inf_aver.columns = adata_sc.uns['mod']['factor_names']

intersect = np.intersect1d(st_exp.var_names, inf_aver.index)
st_exp = st_exp[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()
# prepare anndata for cell2location model
scvi.data.setup_anndata(adata=st_exp)
scvi.data.view_anndata_setup(st_exp)

# create and train the model
mod = cell2location.models.Cell2location(
    st_exp, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=30,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection (using default here):
    detection_alpha=200
)
mod.train(max_epochs=30000,
          # train using full data (batch_size=None)
          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1,
          use_gpu=True)

st_exp = mod.export_posterior(
    st_exp, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
)

# Save model
mod.save(f"{run_name}", overwrite=True)
st_exp.obsm['means_cell_abundance_w_sf'].to_csv(f'{results_folder}/mean.decon.tsv',sep = '\t',header=True,index=True)
st_exp.obsm['stds_cell_abundance_w_sf'].to_csv(f'{results_folder}/stds.decon.tsv',sep = '\t',header=True,index=True)
st_exp.obsm['q05_cell_abundance_w_sf'].to_csv(f'{results_folder}/q05.decon.tsv',sep = '\t',header=True,index=True)
st_exp.obsm['q95_cell_abundance_w_sf'].to_csv(f'{results_folder}/q95.decon.tsv',sep = '\t',header=True,index=True)
# mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)
# Save anndata object with results
# adata_file = f"{run_name}/sp.h5ad"
# st_exp.write(adata_file)
