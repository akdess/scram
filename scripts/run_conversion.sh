#!/bin/bash

# Usage: ./run_conversion.sh <seurat_obj_path>  <output_dir>

SEURATOBJ=$1
OUTPUT_DIR=$2

# Run the R script
Rscript convert_to_h5ad.R $SEURATOBJ $OUTPUT_DIR

# Run the Python script
python3 convert_to_h5ad.py $OUTPUT_DIR 
