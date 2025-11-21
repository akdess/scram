# convert_to_h5ad.R
args <- commandArgs(trailingOnly = TRUE)

# Arguments
seuratObj_path <- args[1]
output_dir <- args[2]

# Load libraries
library(Seurat)
library(Matrix)

# Load data
load(seuratObj_path) # Modify according to your data format

# Create the directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
seuratObj@meta.data$cell_id <- rownames(seuratObj@meta.data)

# Set the working directory to the output directory
setwd(output_dir)

# Create Seurat object
# Write metadata to CSV
write.csv(seuratObj@meta.data, file=paste0(output_dir, "/metadata.csv"), quote=FALSE, row.names=FALSE)

# Get counts matrix and write to MatrixMarket format
counts_matrix <- GetAssayData(seuratObj, assay='RNA', layer='counts')
writeMM(counts_matrix, file=paste0(output_dir, "/counts.mtx"))

# Write gene names to CSV
write.table(data.frame('gene' = rownames(counts_matrix)), file=paste0(output_dir, "/gene_names.csv"), quote=FALSE, row.names=FALSE, col.names=FALSE)
