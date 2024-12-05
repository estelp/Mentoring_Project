#!/bin/bash

############ SLURM Configuration ##############

### Define Job name
#SBATCH --job-name=snp2_calling

### Define partition to use
#SBATCH -p normal

### Define number of CPUs to use
#SBATCH -c 16

### Specify the node to run on
#SBATCH --nodelist=node20  # Run the job on node20

#################################################

########### Execution Commands ###################

# Variables
SORTED_PATH="/scratch/MOryzae/MAPPING/bam_mapped_sort"
VCF_OUTPUT="/scratch/MOryzae/SNP2/vcf_files/all_samples.vcf.gz"
SNP_OUTPUT="/scratch/MOryzae/SNP2/vcf_files/snp.vcf.gz"
REF_GENOME="/scratch/MOryzae/REF/MOryzae_genomic.fna"
SNP_STATS_DIR="/scratch/MOryzae/SNP2/stats"

# Load necessary modules
module load samtools/1.18
module load bcftools/1.18

# Create necessary directories
mkdir -p /scratch/MOryzae/SNP2/vcf_files "$SNP_STATS_DIR"

# Step 1: Generate VCF file, compress, and index
echo -e "######################\nGenerating compressed VCF file"
bcftools mpileup -Ou --threads 16 -f "$REF_GENOME" "$SORTED_PATH"/*.mappedpaired.sorted.bam | \
  bcftools call -mv -Oz -o "$VCF_OUTPUT" || {
    echo "Error: VCF generation and compression failed" >&2
    exit 1
}

# Step 2: Index the compressed VCF file
echo -e "######################\nIndexing compressed VCF file"
bcftools index "$VCF_OUTPUT" || {
    echo "Error: VCF file indexing failed" >&2
    exit 1
}

# Step 3: Filter to retain only SNPs
echo -e "######################\nFiltering to retain only SNPs"
bcftools view -v snps -Oz -o "$SNP_OUTPUT" "$VCF_OUTPUT" || {
    echo "Error: SNP filtering failed" >&2
    exit 1
}

# Step 4: Index the SNP VCF file
echo -e "######################\nIndexing SNP VCF file"
bcftools index "$SNP_OUTPUT" || {
    echo "Error: SNP VCF file indexing failed" >&2
    exit 1
}

# Step 5: Generate SNP statistics
echo -e "######################\nGenerating SNP statistics"
bcftools stats "$SNP_OUTPUT" > "$SNP_STATS_DIR/all_samples_SNP_statistics.txt" || {
    echo "Error: Failed to generate SNP statistics" >&2
    exit 1
}

echo "Compressed VCF file and SNP-specific VCF file generated and indexed successfully."
