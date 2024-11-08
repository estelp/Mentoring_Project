#!/bin/bash

############# SLURM Configuration ##############

### Define Job name
#SBATCH --job-name=multiqc

### Define partition to use
#SBATCH -p normal

### Define number of CPUs to use
#SBATCH -c 8

### Specify the node to run on
#SBATCH --nodelist=node13  # Specifies that the job should run on node13

#################################################

########### Execution Command ##################

# Define paths to the working directories
FASTQC_PATH="/scratch/MOryzae/QC/FastQC"
MULTIQC_OUTPUT_PATH="/scratch/MOryzae/QC/MultiQC"

# Load the MultiQC module
module load multiqc/1.9

# Create the output directory
mkdir -p "$MULTIQC_OUTPUT_PATH"

# Run MultiQC on the FastQC reports
echo "Running MultiQC on FastQC reports in $FASTQC_PATH..."
multiqc "$FASTQC_PATH" -o "$MULTIQC_OUTPUT_PATH"

# Check if MultiQC ran successfully
if [[ $? -eq 0 ]]; then
    echo "MultiQC completed successfully."
else
    echo "Error: MultiQC encountered an issue."
fi
