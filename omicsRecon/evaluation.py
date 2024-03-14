import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def cal_cell_cor(filter_pred_rna, filter_true_rna, shared_idx):
    '''
    This function estimate the gene expression correlation between two df.
    '''
    cor_lst = []
    for i in range(filter_pred_rna.shape[0]):
        cor = np.corrcoef(filter_pred_rna.iloc[i],filter_true_rna.iloc[i])[0,1]
        cor_lst.append(cor)
    tmp = pd.DataFrame(cor_lst,columns = ['cor'],index = shared_idx)
    return tmp


def cal_gene_cor(filter_pred_rna, filter_true_rna, shared_gene):
    '''
    This function estimate the gene expression correlation between two df.
    '''
    cor_lst = []
    for gene in filter_pred_rna.columns:
        cor = np.corrcoef(filter_pred_rna[gene],filter_true_rna[gene])[0,1]
        cor_lst.append(cor)
    tmp = pd.DataFrame(cor_lst,columns = ['cor'],index = shared_gene)
    return tmp


def cell_cor(pred_rna,true_rna):
    shared_idx = pred_rna.index.intersection(true_rna.index)
    shared_gene = pred_rna.columns.intersection(true_rna.columns)
    filter_pred_rna = pred_rna.loc[shared_idx, shared_gene]
    filter_true_rna = true_rna.loc[shared_idx, shared_gene]

    cor_df = cal_cell_cor(filter_pred_rna, filter_true_rna, shared_idx)
    
    return cor_df


def gene_cor(pred_rna,true_rna):
    shared_idx = pred_rna.index.intersection(true_rna.index)
    shared_gene = pred_rna.columns.intersection(true_rna.columns)
    filter_pred_rna = pred_rna.loc[shared_idx, shared_gene]
    filter_true_rna = true_rna.loc[shared_idx, shared_gene]

    cor_df = cal_gene_cor(filter_pred_rna, filter_true_rna, shared_gene)
    
    return cor_df
