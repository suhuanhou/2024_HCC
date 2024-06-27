# conda activate scVelo
# python
################################################################################
##### setting(for python)
import anndata
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import re
from pathlib import Path

dir_main = Path("/share/home/shh/HCC_snRNA-seq")
dir_in = dir_main.joinpath("5.Hepatocyte/dataset")
dir_dataset = dir_main.joinpath("26.scVelo/dataset")
dir_result = dir_main.joinpath("26.scVelo/result")
os.chdir(str(dir_result))
################################################################################
##### read data
sample_obs = pd.read_csv("cellID_obs.csv")
umap = pd.read_csv("cell_embeddings.csv")
cell_clusters = pd.read_csv("clusters_obs.csv")

# Filter out Cell IDs 
adata  = anndata.read_loom(dir_dataset.joinpath("merge.loom"))
adata.obs.index = adata.obs.index.map(lambda x: re.sub(r'cellsnp_.*:', '', x))
adata.obs.index = adata.obs.index.str.replace('x', '-1')
adata = adata[np.isin(adata.obs.index, sample_obs["x"])]


# Add umap coordinates and celltype to the anndata object
adata_index = pd.DataFrame(adata.obs.index)
adata_index = adata_index.rename(columns = {"CellID":'Cell ID'})


# cell_clusters=cell_clusters.iloc[:,1:]
cell_clusters = cell_clusters.rename(columns = {'Unnamed: 0':'Cell ID'})
cell_clusters = cell_clusters[np.isin(cell_clusters["Cell ID"],adata_index["Cell ID"])]
cell_clusters_ordered = adata_index.merge(cell_clusters, on = "Cell ID")
adata.obs['clusters'] = cell_clusters_ordered.set_index('Cell ID')['x']


bdata = adata.copy() # adata = bdata.copy()
umap = umap.rename(columns = {'Unnamed: 0':'Cell ID'})
umap = umap[np.isin(umap["Cell ID"], adata_index["Cell ID"])] 
umap_ordered = adata_index.merge(umap, on = "Cell ID")
umap_ordered = umap_ordered.iloc[:,1:]
adata.obsm['X_umap'] = umap_ordered.values 


# Data preprocessing
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)
scv.tl.velocity(adata, mode = "stochastic")
scv.tl.velocity_graph(adata)


##### visualization
scv.pl.velocity_embedding_stream(adata, basis="umap", legend_loc='none', dpi=600, save="umap.png", palette = ['#3EA94B','#3E86BD','#E51F17','#DE75DD','#EC8392','#31AADF','#FDCE4F','#EC7B1A','#1505A8','#C20576'])

# spliced/unspliced
scv.pl.proportions(adata, save = True)