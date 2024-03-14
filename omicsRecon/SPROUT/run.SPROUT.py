import sys
sys.path.insert(1, '/public/wanwang6/5.Simpute/SPROUT_fast/')
import os
import scanpy as sc
import pandas as pd
import numpy as np
import json
from src import sprout
import warnings
warnings.filterwarnings("ignore")
lr_dir = '/public/wanwang6/5.Simpute/SPROUT_fast/LR/'
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-s', '--st-file', dest='st_file', required=True, help='Spatial transcriptomics data')
parser.add_argument('-c', '--st-coordinate', dest='st_coord', required=False, help='Spatial coordinates of the spatial transcriptomics data')
parser.add_argument('-v', '--weight-file', dest='w_file', required=False, help='Deconvoluted ST data or software (tangram or c2l). If not provided, use c2l.')
parser.add_argument('-r', '--sc-file', dest='sc_file', required=True, help='Single-cell candidate library of the corresponding ST tissue')
parser.add_argument('-m', '--meta-file', dest='meta_file', required=True, help='Cell-type annotation of the single-cell candidate library')
parser.add_argument('-a', '--is-human', dest='is_human', required=True, default=True,help='If the species is human, default True')
parser.add_argument('-t', '--type-key', dest='tp_key', required=True, help='The colname of celltype in SC_meta')
parser.add_argument('-b', '--batch-key', dest='batch_key', required=False, help='batch column name in SC meta')
parser.add_argument('-p','--cell-num-per-spot', dest='num_per_spot', required=False,default=10, help='Estimated cell number per spot. Default is 10')
parser.add_argument('-n', '--name', dest='name', required=False, help='Sample name which will be set as the prefix of output')
parser.add_argument('-o', '--out-dir', dest='out_dir', required=False, help='Output file path')
args = parser.parse_args()

def read_csv_tsv(filename):
    if 'csv' in filename:
        tmp = pd.read_csv(filename, sep = ',',header = 0,index_col=0)
    else:
        tmp = pd.read_csv(filename, sep = '\t',header = 0,index_col=0)
    return tmp
# st_decon = read_csv_tsv('/public/wanwang6/5.Simpute/7.m_embryo/2.results/1.tangram/tg_st_weight.tsv')
# st_exp = read_csv_tsv('/public/wanwang6/5.Simpute/7.m_embryo/1.input/1.E15/E15_ST_1853.tsv')
# sp = 'mouse'
# lr_df = pd.read_csv(f'{lr_dir}/{sp}_LR_pairs.txt',sep='\t',header=None)
# sp = 'human' if is_human else 'mouse'
st_exp = read_csv_tsv(args.st_file)
st_coord = read_csv_tsv(args.st_coord)
sc_exp = read_csv_tsv(args.sc_file)
sc_meta = read_csv_tsv(args.meta_file)
sc_meta['celltype'] = sc_meta[args.tp_key]
sc_exp = sc_exp.loc[sc_meta.index]
sp = args.is_human
lr_df = pd.read_csv(f'{lr_dir}/{sp}_LR_pairs.txt',sep='\t',header=None)
print(f'Using {sp} LRdb...')
# 1. outdir
if args.out_dir is not None:
    outDir = args.out_dir + '/' 
    if not os.path.exists(outDir):
        os.makedirs(outDir)
else:
    outDir = os.getcwd() + '/'
# 2.sample name as prefix
if args.name is not None:
    outDir = f'{outDir}/{args.name}'
    if outDir[-1] == '/':
        # name end with / creat new dir
        if not os.path.exists(outDir):
            os.makedirs(outDir)

# 4.st meta
if 'row' in st_coord.columns.tolist():
    # if the x and y coordinates are indexed with row and col
    st_coord = st_coord[['row','col']].copy()
else:
    st_coord = st_coord[['x','y']].copy()
st_coord.columns = ['x','y']
######################### decon software #########################
if (args.w_file is None) or (args.w_file == 'c2l'):
    # if haven't run c2l yet
    if not os.path.exists(f'{outDir}/c2l/q05.decon.tsv'):
        if not os.path.exists(f'{outDir}/c2l/'):
            os.makedirs(f'{outDir}/c2l/')
        print(f'Saving the c2l results in {outDir}/c2l/')
        # run c2l
        import subprocess
        print('Running c2l...')
        subprocess.run(["/home/grads/wanwang6/anaconda3/envs/cell2loc_env/bin/python", "/public/wanwang6/5.Simpute/4.compare_softwares/a.SPROUT/0.scripts/c2l.py", \
            '-s',args.st_file, \
            '-c',args.st_coord, \
            '-r',args.sc_file, \
            '-m',args.meta_file, \
            '-t',args.tp_key, \
            '-b',args.batch_key, \
            '-o',f'{outDir}/c2l/'])
    st_decon = pd.read_csv(f'{outDir}/c2l/q05.decon.tsv',header=0, index_col=0, sep = '\t')
    st_decon.columns = [x.split('sf_')[-1] for x in st_decon.columns]
    st_decon = pd.DataFrame(st_decon).div(np.sum(st_decon, axis = 1), axis = 0)
elif args.w_file == 'tg':
    # if haven't run tangram yet
    if not os.path.exists(f'{outDir}/tg/magic_tg_st_weight.tsv'):
        if not os.path.exists(f'{outDir}/tg/'):
            os.makedirs(f'{outDir}/tg/')
        print(f'Saving the tangram results in {outDir}/tg/')
        # run tangram
        import subprocess
        print('Running tangram...')
        subprocess.run(["/home/grads/wanwang6/anaconda3/envs/tangram-env/bin/python", \
        "/public/wanwang6/5.Simpute/4.compare_softwares/4.Tangram/0.script/0.tangram.py", \
            '-s',args.st_file, \
            '-c',args.st_coord, \
            '-r',args.sc_file, \
            '-m',args.meta_file, \
            '-t',args.tp_key, \
            '-a',args.is_human, \
            '-o',f'{outDir}/tg/'])
    st_decon = pd.read_csv(f'{outDir}/tg/tg_st_weight.tsv',header=0, index_col=0, sep = '\t')
else:
    if 'c2l' in args.w_file:
        st_decon = read_csv_tsv(args.w_file)
        st_decon.columns = [x.split('sf_')[-1] for x in st_decon.columns]
        st_decon = pd.DataFrame(st_decon).div(np.sum(st_decon, axis = 1), axis = 0)
    else:
        st_decon = read_csv_tsv(args.w_file)
######################### decon software #########################
st_decon = pd.DataFrame(st_decon).div(np.sum(st_decon, axis = 1), axis = 0)
st_exp.index = st_exp.index.map(str)
st_decon.index = st_decon.index.map(str)
st_coord.index = st_coord.index.map(str)
num_per_spot = int(args.num_per_spot)
repeat_penalty = np.round((st_exp.shape[0] * num_per_spot/sc_exp.shape[0]) * 10)
print(f'Setting {num_per_spot} cells per spot; half when chosen {repeat_penalty} times...')
sprout_obj = sprout.SPROUT(st_exp = st_exp, st_coord = st_coord, weight = st_decon, 
        sc_exp = sc_exp, meta_df = sc_meta, cell_type_key = 'celltype',lr_df = lr_df,
        save_path = outDir)
spot_cor,picked_index_df = sprout_obj.select_sc(num_per_spot = num_per_spot, mode = 'strict', max_rep = 1, repeat_penalty = repeat_penalty)
sc_coord = sprout_obj.spatial_recon(left_range = 0, right_range = 10, steps = 1, dim = 2,max_dist = 1)
sc_agg_exp = sc_exp.loc[sprout_obj.picked_index_df['sc_id']]
sc_agg_exp = sc_agg_exp.reset_index()
# delet original index
del sc_agg_exp[sc_agg_exp.columns[0]]
sc_agg_exp.to_csv(f'{outDir}/sc_agg_exp.tsv',sep = '\t',header=True,index=True)
