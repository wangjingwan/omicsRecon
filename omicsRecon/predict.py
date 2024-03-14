import os
import pandas as pd

def atac2RNA(pred_method = 'signac', species = None, Rscript_path = None, outDir = None,
            ATAC_h5Seurat = None, frag_Seurat = None, 
            ATAC_rds = None):
    '''
            @author: Jingwan WANG
            This function estimate the gene expression from single-cell ATAC data.
                        
            Parameters
            ----------
            pred_method : str
            prediction method, currently only support seurat.

            species : str
            species of the data, currently only support human and mouse.

            ATAC_h5Seurat : str
            ATAC file with h5Seurat format. No obj@assays$ATAC@fragments yet.
            
            frag_Seurat : str
            fragment file for Seurat. Normally names atac_fragments.tsv.gz
            
            ATAC_rds : str
            ATAC file with rds format. Comparing with ATAC_h5Seurat, normally already have obj@assays$ATAC@fragments in this seurat object.

            Rscript_path : str
            which Rscript path used for seurat

            outDir : str
            Output path
    '''
    if pred_method == 'signac':
        if not os.path.exists(f'{outDir}/atac_pred_RNA.tsv'):
            if ATAC_h5Seurat:
                if frag_Seurat:
                    # both provided
                    script_path = os.path.dirname(os.path.realpath(__file__)) + '/scripts/ATAC_seurat_pipeline.R'
                    os.system(str(f'{Rscript_path} --vanilla {script_path} -s {species} -a {ATAC_h5Seurat} -f {frag_Seurat} -o {outDir}'))
                else:
                    raise ValueError('Please provide the paired fragment file for Seurat.')

            if ATAC_rds:
                script_path = os.path.dirname(os.path.realpath(__file__)) + '/scripts/ATAC_seurat_pipeline.R'
                os.system(str(f'{Rscript_path} --vanilla {script_path} -s {species} --rds_path {ATAC_rds} -o {outDir}'))
                # os.system(str(f'{Rscript_path} --vanilla /data6/wangjingwan/0.omics/sprout_omics-main/omicsRecon/scripts/ATAC_seurat_pipeline.R -s {species} --rds_path {ATAC_rds} -o {outDir}'))
            # if both none
            if not ATAC_h5Seurat and not ATAC_rds:
                raise ValueError('Please provide the ATAC file with h5Seurat format or rds format.')
            
    if pred_method == 'archr':
        print('Building...')
    # load pred_rna
    pred_rna = pd.read_csv(f'{outDir}/atac_pred_RNA.tsv',header=0,index_col=0,sep='\t')
    return pred_rna
    




# class OmicsToRNA:
#     def __init__(self, RNA_file = None, RNA_h5ad = None, RNA_h5Seurat = None, 
#                  ATAC_h5Seurat = None, frag_Seurat = None, species = None,
#                  train_methy = None, train_rna = None, test_methy = None, test_rna = None,
#                  celltype_key = None, celltype_file = None, 
#                  python_path = None, Rscript_path = None, outDir = None):
#         """
#             @author: Jingwan WANG
#             This function estimate the gene expression from single-cell omics data.
                        
#             Parameters
#             ----------
            
#             RNA_file : str
#             scRNA-seq data count file.
            
#             RNA_h5ad : str
#             scRNA-seq data file with h5ad format.
            
#             RNA_h5Seurat : str
#             scRNA-seq data file with h5Seurat format.
                        
#             ATAC_h5Seurat : str
#             ATAC file with h5Seurat format.
            
#             frag_Seurat : str
#             fragment file for Seurat. Normally names atac_fragments.tsv.gz

#             train_methy : str
#             scBS-seq data beta score file for training.

#             train_rna : str
#             Matching scRNA-seq data count file for training.

#             test_methy : str
#             scBS-seq data beta score file for prediction.

#             test_rna : str
#             scRNA-seq data count file for prediction evaluation [optional].

#             celltype_key : str
#             celltype annotataion title in scRNA-seq data h5ad file or h5Seurat file
            
#             celltype_file : str
#             celltype annotataion file
            
#             python_path : str
#             which python path used for DL
            
#             Rscript_path : str
#             which Rscript path used for seurat

#             outDir : str
#             Output path
            
#             """
        
#         self.RNA_file = RNA_file
#         self.RNA_h5ad = RNA_h5ad
#         self.RNA_h5Seurat = RNA_h5Seurat
#         self.ATAC_h5Seurat = ATAC_h5Seurat 
#         self.frag_Seurat = frag_Seurat
#         self.species = species
#         self.train_methy = train_methy
#         self.train_rna = train_rna
#         self.test_methy = test_methy
#         self.test_rna = test_rna

#         self.celltype_key = celltype_key
#         self.celltype_file = celltype_file
#         self.python_path = python_path
#         self.Rscript_path = Rscript_path
#         self.outDir = outDir
    
#     def transform(self, dataform, method):
#         if self.python_path is None:
#             self.python_path = 'python'
#         if self.Rscript_path is None:
#             self.Rscript_path = 'Rscript'
#         # ATAC
#         if dataform == 'atac':
#             if method == 'seurat':
#                 h5_path = self.h5_path
#                 frag_path = self.frag_path
#                 outDir = self.outDir
#                 species = self.species
#                 Rscript_path = self.Rscript_path
#                 os.system(f'{Rscript_path} ../calledScripts/ATAC_seurat_pipeline.R {h5_path} {frag_path} {outDir} {species}')

#         # methy
#         if dataform == 'methy':
#             if method =='dl':
#                 train_methy = self.train_methy
#                 train_rna = self.train_rna
#                 test_methy = self.test_methy
#                 test_rna = self.test_rna
#                 celltype_file = self.celltype_file
#                 celltype_key = self.celltype_key
#                 outDir = self.outDir
#                 species = self.species
#                 python_path = self.python_path
#                 os.system(f'{python_path} ../calledScripts/methy_pipeline.py {train_methy} {train_rna} {test_methy} {test_rna} {celltype_file} {celltype_key} {outDir}')
#         # protein
#         if dataform =='protein':
#             h5_path = self.h5_path
#             frag_path = self.frag_path
#             outDir = self.outDir
#             species = self.species
#             Rscript_path = self.Rscript_path
#             os.system(f'{Rscript_path} ../calledScripts/ATAC_seurat_pipeline.R {h5_path} {frag_path} {outDir} {species}')
