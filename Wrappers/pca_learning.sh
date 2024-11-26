#!/bin/bash

############# SLURM Configuration ##############
#SBATCH --job-name=pca_learning
#SBATCH -p normal
#SBATCH -c 8
#SBATCH --nodelist=node20

#################################################

# Load necessary modules
module load python/3.12.0

# Define directories
PCA_RESULTS_DIR="/scratch/MOryzae/PLINK"
OUTPUT_PLOT_DIR="/scratch/MOryzae/PCA"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_PLOT_DIR"

# List of subdirectories to process
DIRECTORIES=("all_samples" "only_snp")

# Loop over each directory
for DIRECTORY in "${DIRECTORIES[@]}"; do
    PCA_FILE="$PCA_RESULTS_DIR/$DIRECTORY/dataset.eigenvec"
    OUTPUT_DIR="$OUTPUT_PLOT_DIR/$DIRECTORY"

    # Check if the PCA file exists
    if [[ -f "$PCA_FILE" ]]; then
        mkdir -p "$OUTPUT_DIR"
        echo "Processing PCA file: $PCA_FILE"

        # Run the Python script for clustering and plotting
        python3 clustering_analysis.py "$PCA_FILE" "$OUTPUT_DIR"

        echo "Processing completed for directory: $DIRECTORY"
    else
        echo "PCA file not found: $PCA_FILE"
    fi
done

