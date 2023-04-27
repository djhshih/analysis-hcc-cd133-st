#!/usr/bin/env python3

import scanpy as sc
import numpy as np
import matplotlib.pyplot as pl
import matplotlib as mp
import cell2location as cl
import seaborn as sb
import pandas as pd
import sys

mp.rcParams['pdf.fonttype'] = 42

outdir = 'immune'

# sample = 'DTA17'
# sample = 'DTA23'
# sample = 'DTA18'
# sample = 'DTA24'

sample = sys.argv[1]
print(sample)

outdir_ref_sig = f'{outdir}/ref_signatures'
outdir_map = f'{outdir}/cell2location_map/{sample}'

# reference signatures
ref_sig = pd.read_parquet(f'{outdir}/ref_signatures/ref_sig.pq')

# input visium data
adata_vis = sc.read_visium(f'../spaceranger/{sample}/outs')

# manually annotated regions
regions = pd.read_csv(
    f'../spaceranger/{sample}/outs/{sample}_region.csv',
    index_col = 'Barcode',
)

# set up sample name
adata_vis.obs['sample'] = list(adata_vis.uns['spatial'].keys())[0]

# add region annotations
adata_vis.obs = adata_vis.obs.join(regions, how='left')

# remove duplicate gene names
dup = adata_vis.var_names.duplicated()
adata_vis = adata_vis[:, ~dup].copy()

# find mitochondria genes
is_mito = np.array([
    gene.startswith('mt-')
    for gene in adata_vis.var_names
])

# backup mitochondria gene expression
adata_vis.obsm['MT'] = adata_vis[:, is_mito].X.toarray()

# remove mitochondria genes
adata_vis = adata_vis[:, ~is_mito]

# subset common genes

genes = np.intersect1d(adata_vis.var_names, ref_sig.index)
len(genes)
pd.Series(genes).duplicated().sum()

adata_vis = adata_vis[:, genes].copy()
ref_sig = ref_sig.loc[genes, :].copy()

# model fitting
cl.models.Cell2location.setup_anndata(
    adata = adata_vis,
    batch_key = 'sample'
)

mod_cl = cl.models.Cell2location(
    adata_vis,
    cell_state_df = ref_sig,
    N_cells_per_location = 10,
    detection_alpha = 20
)
mod_cl.view_anndata_setup()

mod_cl.train(
    max_epochs = 30000,
    batch_size = None,   # use full data
    train_size = 1,
    use_gpu = True,
)

adata_vis = mod_cl.export_posterior(
    adata_vis,
    sample_kwargs = {
        'num_samples': 1000,
        'batch_size': mod_cl.adata.n_obs,
        'use_gpu': True
    }
)

pl.axline((0, 0), slope=1, color='grey')
mod_cl.plot_QC()
pl.show()

# add 5% quantile
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = \
    adata_vis.obsm['q05_cell_abundance_w_sf']

# save model
mod_cl.save(outdir_map, overwrite=True)
adata_vis.write(f'{outdir_map}/{sample}.h5ad')

