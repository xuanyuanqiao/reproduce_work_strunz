## figure 6a (It was originally run in jupyter)
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import scvelo as scv

adata = sc.read("/mnt/c/2019/data/EpiHiRes.h5ad",
                cache = True)

sc.settings.verbosity = 3            
sc.logging.print_header()
sc.settings.set_figure_params(dpi=150, facecolor='white')

##plot the velocity
scv.pl.velocity_embedding_grid(adata, basis = "umap",
                               color = "louvain",title='RNA velocity',
                               scale=0.3,alpha=0.05,save=True)

##figure 4b
sc.pl.umap(adata, color="cell_type_recombined",
           title='cell type',edges=False,save=True)

##figure 4c
sc.pl.umap(adata, color = "grouping",
           save=True,title='time after bleo')

##figure 4e
##calculate PAGA
sc.tl.paga(adata, groups = "louvain")
sc.pl.paga(adata, threshold=0.3, edge_width_scale=.2,save=True)

##use PAGA compare
sc.pl.paga_compare(
    adata, threshold=0.3, title='PAGA', right_margin=0.2,
    size=10, edge_width_scale=0.2,
    legend_fontsize=5, fontsize=8, frameon=False, edges=True,save=True)
