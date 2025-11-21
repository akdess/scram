# convert_to_h5ad.py
import sys
import scipy.io
import scanpy as sc
import pandas as pd

# Arguments
input_dir = sys.argv[1]

# File paths
counts_path = f'{input_dir}/counts.mtx'
metadata_path = f'{input_dir}/metadata.csv'
gene_names_path = f'{input_dir}/gene_names.csv'

# Read the MatrixMarket file in Python
X = scipy.io.mmread(counts_path)

# Create AnnData object
adata = sc.AnnData(X=X.transpose().tocsr())

# Load cell metadata
cell_meta = pd.read_csv(metadata_path)

# Load gene names
with open(gene_names_path, 'r') as f:
    gene_names = f.read().splitlines()

# Assign metadata and gene names to AnnData
adata.obs = cell_meta
adata.var.index = gene_names
adata.obs_names = adata.obs["cell_id"] 

output_file = f'{input_dir}/adata.h5ad'

# Save the AnnData object as h5ad
adata.write(output_file)


