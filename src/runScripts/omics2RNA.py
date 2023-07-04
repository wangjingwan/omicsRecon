import sys
import os

class OmicsToRNA:
    def __init__(self, RNA_file = None, RNA_h5ad = None, RNA_h5Seurat = None, Spatial_file = None, Spatial_h5ad = None, Spatial_h5Seurat = None, celltype_key = None, celltype_file = None, python_path = None, output_path = None):
        """
            @author: Jingwan WANG
            This function integrates spatial and scRNA-seq data to predictes the celltype deconvolution of the spots.
            
            A minimal example usage:
            Assume we have (1) scRNA-seq data file named RNA_h5ad or RNA_h5Seurat
            (2) spatial transcriptomics data file named Spatial_h5ad or Spatial_h5Seurat
            (3) celltype annotataion title in scRNA-seq data file
            
            >>> import Benchmarking.DeconvolutionSpot as DeconvolutionSpot
            >>> test = DeconvolutionSpot.Deconvolutions(RNA_file, RNA_h5ad, RNA_h5Seurat, Spatial_file, Spatial_h5ad, Spatial_h5Seurat, celltype_key, celltype_file, python_path, output_path)
            >>> Methods = ['Cell2location','SpatialDWLS','RCTD','STRIDE','Stereoscope','Tangram','DestVI', 'Seurat', 'SPOTlight', 'DSTG']
            >>> Result = test.Dencon(Methods)
            
            Parameters
            -------
            
            RNA_file : str
            scRNA-seq data count file.
            
            RNA_h5ad : str
            scRNA-seq data file with h5ad format.
            
            RNA_h5Seurat : str
            scRNA-seq data file with h5Seurat format.
            
            Spatial_file : str
            Spatial data count file.
            
            Spatial_h5ad : str
            Spatial data file with h5ad format.
            
            Spatial_h5Seurat : str
            Spatial data file with h5Seurat format.
            
            celltype_key : str
            celltype annotataion title in scRNA-seq data h5ad file or h5Seurat file
            
            celltype_file : str
            celltype annotataion file
            
            python_path : str
            which python path used for Cell2location
            
            output_path : str
            Outfile path
            
            """
        
        self.RNA_file = RNA_file
        self.RNA_h5ad = RNA_h5ad
        self.RNA_h5Seurat = RNA_h5Seurat
        self.Spatial_file = Spatial_file
        self.Spatial_h5ad = Spatial_h5ad
        self.Spatial_h5Seurat = Spatial_h5Seurat
        self.celltype_key = celltype_key
        self.celltype_file = celltype_file
        self.python_path = python_path
        self.output_path = output_path
    
    def transform(self, dataform, method):
        if dataform =='protein':
            RNA_h5Seurat = self.RNA_h5Seurat
            Spatial_h5Seurat = self.Spatial_h5Seurat
            celltype_key = self.celltype_key
            output_path = self.output_path
            python_path = self.python_path
            os.system('Rscript ../calledScripts/ATAC_seurat_pipeline.R ' + RNA_h5Seurat + ' ' + Spatial_h5Seurat + ' ' + celltype_key + ' ' +  python_path + ' ' + output_path)

        if dataform == 'atac':
            if method == 'seurat':
                RNA_h5ad = self.RNA_h5ad
                Spatial_h5ad = self.Spatial_h5ad
                celltype_key = self.celltype_key
                output_path = self.output_path
                os.system('Rscript ../calledScripts/ATAC_seurat_pipeline.R ' + RNA_h5ad + ' ' + Spatial_h5ad + ' ' + celltype_key + ' ' + output_path)

        if dataform =='methy':
            RNA_h5Seurat = self.RNA_h5Seurat
            Spatial_h5Seurat = self.Spatial_h5Seurat
            celltype_key = self.celltype_key
            output_path = self.output_path
            python_path = self.python_path
            os.system('python Codes/Deconvolution/Cell2location_pipeline.py ' + RNA_h5ad + ' ' + Spatial_h5ad + ' ' + celltype_key + ' ' + output_path)