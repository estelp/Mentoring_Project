#!/bin/bash

############# SLURM Configuration ##############

### Define Job name
#SBATCH --job-name=eigenvec_to_csv

### Define partition to use
#SBATCH -p normal

### Define number of CPUs to use
#SBATCH -c 8

### Specify the node to run on
#SBATCH --nodelist=node20  # Run the job on node20

#################################################

########### Execution Command ###################

# Define directories
INPUT_DIR="/scratch/MOryzae/PLINK2"
OUTPUT_DIR="/scratch/MOryzae/PLINK2"

# List of directories to process
DIRECTORIES=("plink")

# Loop through each directory
for DIR in "${DIRECTORIES[@]}"; do
    PCA_FILE="$INPUT_DIR/$DIR/dataset.eigenvec"
    OUTPUT_CSV="$OUTPUT_DIR/$DIR/dataset.csv"

    # Check if PCA results exist
    if [ -f "$PCA_FILE" ]; then
        echo "Processing PCA results for $DIR..."

        # Convert eigenvec file to CSV format
        awk 'NR==1{print "FID,IID,PC1,PC2,PC3"} NR>1{print $1","$2","$3","$4","$5}' "$PCA_FILE" > "$OUTPUT_CSV"

        # Check if conversion was successful
        if [ $? -eq 0 ]; then
            echo "Conversion successful: $OUTPUT_CSV created."
        else
            echo "Error: Failed to convert $PCA_FILE to CSV."
            exit 1
        fi
    else
        echo "Error: PCA results file not found in $DIR."
        exit 1
    fi
done

echo "PCA analysis completed. Results are saved in $OUTPUT_DIR."

