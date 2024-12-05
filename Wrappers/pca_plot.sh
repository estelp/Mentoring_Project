#!/bin/bash

############# SLURM Configuration ##############

### Define Job name
#SBATCH --job-name=genome_pca_plot

### Define partition to use
#SBATCH -p normal

### Define number of CPUs to use
#SBATCH -c 8

### Specify the node to run on
#SBATCH --nodelist=node20

#################################################

########### Execution Command ###################

module load python/3.12.0  # Charge Python 3.12 sur le cluster

# Define directories
PCA_RESULTS_DIR="/scratch/MOryzae/PLINK"
OUTPUT_PLOT_DIR="/scratch/MOryzae/PLINK"

# List of directories to process
DIRECTORIES=("plink")

# Loop through each directory
for DIRECTORY in "${DIRECTORIES[@]}"; do
    PCA_FILE="$PCA_RESULTS_DIR/$DIRECTORY/dataset.eigenvec"
    OUTPUT_DIR="$OUTPUT_PLOT_DIR/$DIRECTORY"
    OUTPUT_PLOT_2D="$OUTPUT_DIR/dataset_2D.png"
    OUTPUT_PLOT_3D="$OUTPUT_DIR/dataset_3D.png"

    # Ensure the output directory exists
    mkdir -p "$OUTPUT_DIR"

    # Check if the PCA results file exists
    if [ -f "$PCA_FILE" ]; then
        echo "Processing PCA results for $DIRECTORY..."
        
        # Call the Python script to generate the plots
        python3 <<EOF
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

# Define input and output paths
pca_results_file = "$PCA_FILE"
output_plot_2D = "$OUTPUT_PLOT_2D"
output_plot_3D = "$OUTPUT_PLOT_3D"

# Read PCA results
try:
    pca_results = pd.read_csv(pca_results_file, sep=r'\s+', header=None)
    pca_results.columns = ['FID', 'IID', 'PC1', 'PC2', 'PC3']
except Exception as e:
    print(f"Error reading PCA results file {pca_results_file}: {e}")
    exit(1)

# Plot 2D scatter plot for PC1 vs PC2
try:
    plt.figure(figsize=(8, 6))
    plt.scatter(pca_results['PC1'], pca_results['PC2'], s=100)
    plt.title('PCA Results: $DIRECTORY (2D)')
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.grid()
    plt.savefig(output_plot_2D)
    plt.close()
except Exception as e:
    print(f"Error creating 2D plot for {pca_results_file}: {e}")
    exit(1)

# Plot 3D scatter plot for PC1, PC2, and PC3
try:
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    scatter = ax.scatter(pca_results['PC1'], pca_results['PC2'], pca_results['PC3'], s=100, c='blue', alpha=0.7)
    ax.set_title('PCA Results: $DIRECTORY (3D)')
    ax.set_xlabel('Principal Component 1')
    ax.set_ylabel('Principal Component 2')
    ax.set_zlabel('Principal Component 3')
    plt.savefig(output_plot_3D)
    plt.close()
except Exception as e:
    print(f"Error creating 3D plot for {pca_results_file}: {e}")
    exit(1)

# Verify plots were created
if not os.path.exists(output_plot_2D) or not os.path.exists(output_plot_3D):
    print(f"Error: Output plots not created for {pca_results_file}")
    exit(1)
EOF

        # Check if the Python script executed successfully
        if [ $? -eq 0 ]; then
            echo "Plots successfully created for $DIRECTORY: $OUTPUT_PLOT_2D, $OUTPUT_PLOT_3D"
        else
            echo "Error: Failed to create plots for $DIRECTORY."
            exit 1
        fi
    else
        echo "Error: PCA results file not found in $DIRECTORY."
        exit 1
    fi
done

echo "All PCA plots created successfully."

