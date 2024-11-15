#!/bin/bash

############# SLURM Configuration ##############

### Define Job name
#SBATCH --job-name=snp_calling

### Define partition to use
#SBATCH -p normal

### Define number of CPUs to use
#SBATCH -c 16

### Specify the node to run on
#SBATCH --nodelist=node20  # Spécifie que le job doit être exécuté sur node20

#################################################

########### Execution Command ##################

# Define file paths
SORTED_PATH="/scratch/MOryzae/MAPPING/bam_mapped_sort"
BCF_PATH="/scratch/MOryzae/SNP/bcf_files/all_samples.bcf"
VCF_PATH="/scratch/MOryzae/SNP/vcf_files/all_samples.vcf"
REF_GENOME="/scratch/MOryzae/REF/MOryzae_genomic.fna"
SNP_STATS_DIR="/scratch/MOryzae/SNP/stats"
OUTPUT_VCF="${VCF_PATH}.gz"
SNP_FILE="/scratch/MOryzae/SNP/vcf_files/all_samples_snp.vcf"
OUTPUT_VCF2="${SNP_FILE}.gz"
ALLELE_FREQ_PATH="/scratch/MOryzae/SNP/allele_frequence"

# Load necessary modules
module load samtools/1.18
module load bcftools/1.18
module load vcftools/0.1.16

# Create directories if necessary
mkdir -p /scratch/MOryzae/SNP/bcf_files /scratch/MOryzae/SNP/vcf_files "$SNP_STATS_DIR" "$ALLELE_FREQ_PATH"


# Use bcftools for mpileup and variant calling
echo -e "######################\nGenerating BCF file"
bcftools mpileup --threads 16 -f "$REF_GENOME" -O b -o "$BCF_PATH" "$SORTED_PATH"/*.mappedpaired.sorted.bam

echo -e "######################\nVariant calling"
bcftools call --threads 16 -v -c -o "$VCF_PATH" "$BCF_PATH"

echo -e "######################\nGenerating SNP statistics"
bcftools stats "$VCF_PATH" > "$SNP_STATS_DIR/all_samples_SNP_statistics.txt"
# Filter to keep only SNPs
echo -e "######################\nFiltering SNPs"
bcftools view -v snps "$VCF_PATH" -o "$SNP_FILE"

# Compress and index the VCF file
echo -e "######################\nCompressing and indexing VCF file"
bgzip -c "$SNP_FILE" > "$OUTPUT_VCF2"
bcftools index "$OUTPUT_VCF"
bcftools index "$OUTPUT_VCF2"

# Calculate allele frequencies
echo -e "######################\nCalculating allele frequencies"
vcftools --gzvcf "$OUTPUT_VCF" --freq --out "${ALLELE_FREQ_PATH}/AF" --max-alleles 2
vcftools --gzvcf "$OUTPUT_VCF" --freq2 --out "${ALLELE_FREQ_PATH}/AF2" --max-alleles 2

vcftools --gzvcf "$OUTPUT_VCF2" --freq --out "${ALLELE_FREQ_PATH}/AF_2" --max-alleles 2
vcftools --gzvcf "$OUTPUT_VCF2" --freq2 --out "${ALLELE_FREQ_PATH}/AF2_2" --max-alleles 2

echo "BCF and VCF files generated successfully."
