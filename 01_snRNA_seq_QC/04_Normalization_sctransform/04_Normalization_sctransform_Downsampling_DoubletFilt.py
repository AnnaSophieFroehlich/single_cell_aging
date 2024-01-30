# Normalization of single nuclei RNAseq data with sctransform
# run with conda environment sc-norm-sct

# Import packages
import sys
import pandas as pd
import anndata
import scanpy as sc
import numpy as np
import scipy as sp
from scipy.sparse import issparse

# Set pathes
sys.path.insert(0,'..') # add parent directory to system pathes so that sct_norm folder can be found
import paths_downsampling as paths
p = paths.get_paths()
print(p)

# import wrapper/interface functions from sct_norm folder
import sct_norm.sctransform as sct


# 1. Load data
adata = sc.read(p['writedir']+'adata_qc_RNA_downsampling_perCell_filtDoublets.h5ad')

# 2. Run sctransform normalization
# sctransform normalizes data and selects also HVGs
# corrected counts are stored in adata.layers['sct_corrected']
# --> set n_top_genes to 5000
sct.sctransform(adata, n_top_genes = 5000)

# 3. Write data
adata.write(p['writedir']+'adata_normlog_RNA_downsampling_perCell_filtDoublets.h5ad')
