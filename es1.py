import scanpy as sc
import os
import matplotlib.pyplot as plt

import scanpy as sc

# Load the PBMC 3k single-cell RNA-seq dataset from Scanpy's built-in datasets
adata = sc.datasets.pbmc3k()

# Filter out cells with fewer than 200 genes
sc.pp.filter_cells(adata, min_genes=200)

# Filter out genes that are expressed in fewer than 3 cells
sc.pp.filter_genes(adata, min_cells=3)

# Annotate mitochondrial genes for quality control
adata.var["mt"] = adata.var_names.str.startswith("MT-")

# Calculate QC metrics including percent of mitochondrial genes
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

# Filter cells with high gene counts or high mitochondrial content (likely low quality)
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]

# Normalize total counts per cell to 10,000 and log-transform the data
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Identify highly variable genes (HVGs) for downstream analysis
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# Keep only highly variable genes
adata = adata[:, adata.var.highly_variable]

# Scale the data to unit variance and zero mean, clip values to max of 10
sc.pp.scale(adata, max_value=10)

# Perform PCA
sc.tl.pca(adata, svd_solver='arpack')

# Compute the neighborhood graph for the cells
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# Compute the UMAP embedding for visualization
sc.tl.umap(adata)

# Perform clustering using the Leiden algorithm
sc.tl.leiden(adata)

# Plot UMAP colored by cluster (Leiden) and known marker genes
sc.pl.umap(adata, color=["leiden", "CD3D", "MS4A1", "GNLY"])

