#!/bin/bash

#SBATCH -c 16
#SBATCH --nodelist=node20

#################################################

########### Execution Command ##################

# Define file paths
SORTED_PATH="/scratch/MOryzae/MAPPING/bam_mapped_sort"
BCF_PATH="/scratch/MOryzae/SNP/bcf_files/all_samples.bcf"
VCF_PATH="/scratch/MOryzae/SNP/vcf_files/all_samples.vcf"
REF_GENOME="/scratch/MOryzae/REF/MOryzae_genomic.fna"
SNP_STATS_DIR="/scratch/MOryzae/SNP/stats"
SNP_FILE="/scratch/MOryzae/SNP/vcf_files/all_samples_snp.vcf"
OUTPUT_VCF="${VCF_PATH}.gz"
OUTPUT_VCF2="${SNP_FILE}.gz"
ALLELE_FREQ_PATH="/scratch/MOryzae/SNP/allele_frequence"

# Load necessary modules
module load samtools/1.18
module load bcftools/1.18
module load vcftools/0.1.16
module load htslib/1.19

# Create directories if necessary
mkdir -p /scratch/MOryzae/SNP/bcf_files /scratch/MOryzae/SNP/vcf_files "$SNP_STATS_DIR" "$ALLELE_FREQ_PATH"

# Step 1: Generate BCF file
echo -e "######################\nGenerating BCF file"
bcftools mpileup --threads 16 -f "$REF_GENOME" -O b -o "$BCF_PATH" "$SORTED_PATH"/*.mappedpaired.sorted.bam || {
    echo "Error: Failed to generate BCF file" >&2
    exit 1
}

# Step 2: Variant calling
echo -e "######################\nVariant calling"
bcftools call --threads 16 -v -c -o "$VCF_PATH" "$BCF_PATH" || {
    echo "Error: Variant calling failed" >&2
    exit 1
}

# Step 3: Generate SNP statistics
echo -e "######################\nGenerating SNP statistics"
bcftools stats "$VCF_PATH" > "$SNP_STATS_DIR/all_samples_SNP_statistics.txt" || {
    echo "Error: Failed to generate SNP statistics" >&2
    exit 1
}

# Step 4: Filter to keep only SNPs
echo -e "######################\nFiltering SNPs"
bcftools view -v snps "$VCF_PATH" -o "$SNP_FILE" || {
    echo "Error: Filtering SNPs failed" >&2
    exit 1
}

# Step 5: Compress and index the VCF file
echo -e "######################\nCompressing and indexing VCF files"
bgzip -c "$VCF_PATH" > "$OUTPUT_VCF" || {
    echo "Error: Compression of VCF file failed" >&2
    exit 1
}
bgzip -c "$SNP_FILE" > "$OUTPUT_VCF2" || {
    echo "Error: Compression of SNP file failed" >&2
    exit 1
}
bcftools index "$OUTPUT_VCF" || {
    echo "Error: Indexing VCF file failed" >&2
    exit 1
}
bcftools index "$OUTPUT_VCF2" || {
    echo "Error: Indexing SNP file failed" >&2
    exit 1
}

# Step 6: Calculate allele frequencies
echo -e "######################\nCalculating allele frequencies"
vcftools --gzvcf "$OUTPUT_VCF" --freq --out "${ALLELE_FREQ_PATH}/AF" --max-alleles 2 || {
    echo "Error: Calculating allele frequencies failed for $OUTPUT_VCF" >&2
}
vcftools --gzvcf "$OUTPUT_VCF" --freq2 --out "${ALLELE_FREQ_PATH}/AF_2" --max-alleles 2 || {
    echo "Error: Calculating allele frequencies failed for $OUTPUT_VCF" >&2
}

echo "BCF and VCF files generated successfully."
