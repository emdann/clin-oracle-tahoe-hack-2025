from pathlib import Path

import numpy as np
import dask.distributed as dd
import scanpy as sc
import anndata as ad
import h5py
import dask

from collections import Counter
import pandas as pd
from tqdm import tqdm
import dask
import time
from dask.distributed import Client, LocalCluster

import dask

sc.logging.print_header()

if use_gpu:
    import rapids_singlecell as rsc
    SPARSE_CHUNK_SIZE = 100_000
    from dask_cuda import LocalCUDACluster

    import cupy as cp


    preprocessing_gpus="0,1,2,3,4,5,6,7"
    cluster = LocalCUDACluster(CUDA_VISIBLE_DEVICES=preprocessing_gpus,
                                threads_per_worker=25,
                                protocol="tcp")
else:
    SPARSE_CHUNK_SIZE = 100_000
    cluster = LocalCluster(n_workers=16)

client = Client(cluster)
if use_gpu:
    mod = rsc
else:
    mod = sc

from packaging.version import parse as parse_version
import gc

adatas = []
all_highly_variable_genes = []

if parse_version(ad.__version__) < parse_version("0.12.0rc1"):
    from anndata.experimental import read_elem_as_dask as read_dask
else:
    from anndata.experimental import read_elem_lazy as read_dask

# non comprehensive filter
def filter_protein_coding_genes(gene_list):
    """
    Filter a list of gene names to keep only protein-coding genes.
    
    Args:
        gene_list: A list of gene names as strings
        
    Returns:
        A list containing only likely protein-coding gene names
    """
    # Common non-coding gene identifiers/patterns
    non_coding_patterns = [
        'LOC', 'LINC', 'MIR', 'mir-', 'SNOR', 'SNHG', 
        'RNU', 'tRNA', 'rRNA', 'snoR', 'snR', 'lncRNA',
        'pseudo', 'NEAT', 'XIST', 'MALAT', 'HOTAIR'
    ]
    
    protein_coding_genes = []
    
    for gene in gene_list:
        # Skip gene if it contains any of the non-coding patterns
        if not any(pattern in gene.upper() for pattern in non_coding_patterns):
            # Skip ENSG IDs (uncharacterized Ensembl genes)
            if not gene.startswith('ENSG'):
                # Skip mitochondrial genes (MT-*)
                if not gene.upper().startswith('MT-'):
                    protein_coding_genes.append(gene)
            
    return protein_coding_genes

for i in tqdm(range(4)):
    id = str(i + 1)
    PATH = f"data/raw/plate{id}_filt_Vevo_Tahoe100M_WServicesFrom_ParseGigalab.h5ad"

    with h5py.File(PATH, "r") as f:
        adata = ad.AnnData(
            obs=ad.io.read_elem(f["obs"]),
            var=ad.io.read_elem(f["var"]),
        )
        adata.X = read_dask(
            f["X"], chunks=(SPARSE_CHUNK_SIZE, adata.shape[1])
        )

    if use_gpu:
        rsc.get.anndata_to_GPU(adata)
    # 100m filtering
    pass_filter_mask = adata.obs["pass_filter"] == "full"
    adata = adata[pass_filter_mask, :].copy()

    # Filter to keep only protein coding genes before HVG selection
    protein_coding_genes = filter_protein_coding_genes(adata.var_names)
    adata = adata[:, protein_coding_genes].copy()

    mod.pp.normalize_total(adata, target_sum=10_000)
    mod.pp.log1p(adata)
    mod.pp.highly_variable_genes(adata, n_top_genes=8_000)

    highly_variable_genes = set(adata.var_names[adata.var["highly_variable"]])
    all_highly_variable_genes.append(highly_variable_genes)
    adatas.append(adata)

## select the genes appears more than two plates
gene_counts = Counter(gene for genes in all_highly_variable_genes for gene in genes)
selected_genes = {gene for gene, count in gene_counts.items() if count > 2}

with open('selected_genes.txt', 'w') as f:
    for gene in selected_genes:
        f.write(f"{gene}\n")
selected_genes = set([x.strip() for x in open("tahoe/selected_genes.txt")])

adatas_to_merge = []
for i in tqdm(range(4)):
    adata = sc.read_h5ad(f"/home/ubuntu/data/raw/plate{i+1}_filt_Vevo_Tahoe100M_WServicesFrom_ParseGigalab.h5ad")
    # Filter to only include selected genes
    common_genes = list(set(adata.var_names) & selected_genes)
    adata = adata[:, common_genes].copy()
    adatas_to_merge.append(adata)

# Merge all datasets at the end
merged_adata = adatas_to_merge[0].concatenate(adatas_to_merge[1:], join='inner')

# Save the merged dataset
merged_adata.write_h5ad("/home/ubuntu/data/merged.h5ad")
