# Population Genetics and Cryptic Species Assessment in *Magnaporthe oryzae*

[![License: AGPL-3.0](https://img.shields.io/badge/License-AGPL--3.0-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)
[![Pipeline: Reproducible](https://img.shields.io/badge/Pipeline-Reproducible-brightgreen.svg)]()
[![Sequencing: Illumina](https://img.shields.io/badge/Sequencing-Illumina-orange.svg)]()
[![Mapping: BWA-MEM](https://img.shields.io/badge/Mapping-BWA--MEM-orange.svg)]()
[![Variant Calling: GATK](https://img.shields.io/badge/Variant%20Calling-GATK-orange.svg)]()
[![Population Structure: LEA/SNMF](https://img.shields.io/badge/Population%20Structure-LEA%2FSNMF-orange.svg)]()
[![Phylogenetics: RAxML-NG](https://img.shields.io/badge/Phylogenetics-RAxML--NG-orange.svg)]()
[![HPC: SLURM](https://img.shields.io/badge/HPC-SLURM-lightgrey.svg)]()
[![Language: Bash | R](https://img.shields.io/badge/Language-Bash%20%7C%20R-blue.svg)]()
[![Status: Complete](https://img.shields.io/badge/Status-Complete-green.svg)]()

## Overview

This repository documents a comprehensive bioinformatics analysis pipeline for investigating population genetic structure and cryptic species diversity in *Magnaporthe oryzae* (M. Grill) Barr, a globally distributed ascomycete fungus responsible for blast disease across diverse cereal and grass hosts. The project employs whole-genome resequencing data from 89 isolates representing multiple host species to elucidate the relationship between host range and genetic population structure, with implications for understanding pathogen emergence, virulence evolution, and host-specificity mechanisms.

---

## Project Information

### Scope and Duration
- **Project Period:** October 21, 2024 – November 25, 2024
- **Preparatory Training:** September 9, 2024 – October 4, 2024
- **Location:** Abidjan, Côte d'Ivoire
- **Framework:** International Certificate in Bioinformatics and Genomics (CIBiG), Central and West African Virus Epidemiology (WAVE) Initiative

### Project Team

**Mentees:**
- Pakyendou Estel NAME (INERA-WAVE, Burkina Faso) — Primary Analysis & Pipeline Development
- Attolou Raoul AGNIMONHAN (AfricaRice, Ivory Coast) — Collaborative Research Support

**Mentors:**
- Aurore COMTE (Institut de Recherche pour le Développement, IRD)
- Sébastien RAVEL (Centre de Coopération Internationale en Recherche Agronomique pour le Développement, CIRAD)

### Contact Information

For questions regarding methodology, reproducibility, or collaboration:

| Contact | Email | Institution |
|---------|-------|-------------|
| Pakyendou Estel NAME | pakyendou.name@gmail.com | INERA-WAVE, Burkina Faso |
| Attolou Raoul AGNIMONHAN | R.Agnimonhan@cgiar.org | AfricaRice, Ivory Coast |
| Aurore COMTE | aurore.comte@ird.fr | IRD |
| Sébastien RAVEL | sebastien.ravel@cirad.fr | CIRAD |

---

## Scientific Rationale and Objectives

### Background

*Magnaporthe oryzae* is a cosmopolitan pathogen of global agricultural significance, causing blast disease in more than 50 plant species encompassing both cultivated cereals (rice, wheat, barley, millet, oats) and wild grasses (perennial ryegrass, St. Augustine grass). While extensively characterized as a rice pathogen due to its devastating agronomic impact, accumulating evidence suggests that *M. oryzae* populations exhibit substantial genetic differentiation based on host identity, reproduction mode, and fungicide sensitivity.

Population genetic structure in fungal plant pathogens reflects both evolutionary processes (mutation, genetic drift, selection) and epidemiological dynamics (host-driven local adaptation, trade-driven dispersal patterns). Understanding these forces is critical for predicting pathogen emergence risks, designing durable resistance strategies, and managing disease in multi-host agricultural systems.

### Research Questions

This project addresses three primary hypotheses:

1. **Host-Specific Lineage Formation:** Do *M. oryzae* isolates from different host species form reproductively isolated, host-specific lineages, or do they represent panmictic populations capable of crossing host barriers?

2. **Cryptic Species Detection:** Do genetic signatures of reproductive isolation define distinct cryptic species within the traditional *M. oryzae* complex, or is observed genetic variation consistent with within-species population substructure?

3. **Population Genetics-Phenotype Relationships:** Do genetic clusters correlate with known epidemiological traits (e.g., fungicide resistance profiles, spore morphology, virulence patterns)?

### Primary Objectives

- Quantify genome-wide genetic variation among 89 *M. oryzae* isolates spanning multiple host origins
- Infer population genetic structure using model-based clustering and phylogenetic approaches
- Identify genetic markers of population divergence and assess signatures of selection
- Evaluate the strength of reproductive isolation between host-associated populations
- Provide genomic context for future comparative genomics and functional studies

---

## Data

### Sequence Data

**Source:** Central and West African Virus Epidemiology (WAVE) Initiative reference collection

**Sample Composition:**
- Total isolates analyzed: 89
- Sequencing platform: Illumina paired-end (specifics in [data_project_pyri.xlsx](/Topic/Data_project_pyri.xlsx))
- Read depth: Sufficient for population-level SNP discovery
- Host species represented: Rice (*Oryza sativa*), wheat, barley, millet, ryegrass, and St. Augustine grass

**Data Organization:**
Raw sequencing data is organized in the `DATA/Raw_Data` directory with the following structure:

```
DATA/
├── Raw_Data/              # Illumina paired-end FASTQ files
│   ├── *_R1.fastq.gz      # Read 1 (forward)
│   └── *_R2.fastq.gz      # Read 2 (reverse)
└── Contigs/               # Reference-aligned assemblies (subset)
```

### Data Availability

Sequencing data is stored on the Montpellier IRD Cluster (bioinfo-master1.ird.fr) under restricted access. Access inquiries should be directed to the project mentors.

Full results, intermediate datasets, and analysis outputs are provided in the [Results](/Results) directory.

---

## Bioinformatics Pipeline

### Overview

The analysis pipeline implements a standard population genomics workflow for fungal whole-genome resequencing data:

```
Raw Sequencing Data
    ↓
Quality Control (FastQC, MultiQC)
    ↓
Read Mapping (BWA-MEM → Reference Genome)
    ↓
Variant Calling (GATK HaplotypeCaller → VCF)
    ↓
Variant Filtering (Quality & Annotation)
    ↓
Population Genetics Analyses
    ├── Genetic Diversity & Polymorphism Indices
    ├── Population Structure (LEA/SNMFAdmixture)
    ├── Phylogenetic Reconstruction (RAxML/IQtree)
    ├── Linkage Disequilibrium Assessment
    └── Evidence for Selection
    ↓
Statistical Validation & Visualization
```

### Computational Environment

- **HPC System:** Montpellier IRD Cluster (bioinfo-master1.ird.fr)
- **Scheduler:** SLURM (Simple Linux Utility for Resource Management)
- **Primary Language:** Bash (workflow automation), R (statistical analysis)

---

## Pipeline Stages

### 1. Cluster Access and Data Setup

#### Login and Node Allocation

```bash
# Connect to the cluster
ssh username@bioinfo-master1.ird.fr

# Request computational resources (interactive session example)
# Parameters: -p (partition), -c (CPU cores), --pty (pseudo-terminal)
srun -p normal -c 8 --pty bash -i
```

#### Working Directory Initialization

```bash
# Create primary working directory
mkdir -p /scratch/MOryzae
cd /scratch/MOryzae

# Verify location
pwd

# Create subdirectories for data organization
mkdir -p DATA/{Raw_Data,Contigs} \
         QC/FastQC \
         SCRIPTS \
         Results/{Mapping,SNP,Population_Structure,Phylogenetics}
```

#### Data Transfer

```bash
# Copy raw sequencing data from NAS storage
scp -r san:/projects/medium/CIBiG_MOryzae/fastq/*_R* /scratch/MOryzae/DATA/Raw_Data

# Verify transfer completeness
ls -lh /scratch/MOryzae/DATA/Raw_Data/ | wc -l  # Should show 178 files (89 pairs)

# Inspect read quality in one file
zcat /scratch/MOryzae/DATA/Raw_Data/sample_name_R1.fastq.gz | head -20
```

---

### 2. Quality Control

#### 2.1 Per-Read Quality Assessment (FastQC)

**Purpose:** Evaluate Illumina read quality, adapter contamination, and sequence anomalies.

**SLURM Script:** `FastQC.sh`

```bash
#!/bin/bash

############# SLURM Configuration ##############
#SBATCH --job-name=qc_fastq
#SBATCH -p normal
#SBATCH -c 2
#SBATCH --array=0-88%4                    # 89 pairs (0-indexed), max 4 concurrent
#SBATCH --nodelist=node20                 # Optional: specify node
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=4G

#################################################

# Set paths
RAW_DATA_PATH="/scratch/MOryzae/DATA/Raw_Data"
QC_OUTPUT="/scratch/MOryzae/QC/FastQC"
mkdir -p ${QC_OUTPUT}

# Get list of unique samples (remove R1/R2 suffix)
SAMPLES=($(ls ${RAW_DATA_PATH}/*_R1.fastq.gz | xargs basename -a | sed 's/_R1.fastq.gz//'))
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

# Execute FastQC
module load fastqc/0.11.9  # Adjust version as needed
fastqc -q \
       --outdir ${QC_OUTPUT} \
       --threads 2 \
       ${RAW_DATA_PATH}/${SAMPLE}_R1.fastq.gz \
       ${RAW_DATA_PATH}/${SAMPLE}_R2.fastq.gz

echo "FastQC completed for sample: ${SAMPLE}"
```

**Submission:**
```bash
cd /scratch/MOryzae/SCRIPTS
sbatch FastQC.sh
```

#### 2.2 Aggregated Quality Report (MultiQC)

```bash
# After all FastQC jobs complete
module load multiqc
cd /scratch/MOryzae/QC
multiqc --outdir ./MultiQC_Report FastQC/
```

**Output:** `MultiQC_Report/multiqc_report.html` — Interactive summary of quality metrics across all samples.

---

### 3. Read Mapping to Reference Genome

#### Reference Genome

Reference: *M. oryzae* 70-15 strain (NCBI Assembly ID: GCF_000146045.2)

**Pre-processing:**
```bash
# Index the reference genome (one-time operation)
REF_GENOME="/scratch/MOryzae/Reference/MOryzae_70-15.fna"
module load bwa
bwa index ${REF_GENOME}
```

#### Alignment Pipeline (BWA-MEM)

**SLURM Script:** `Mapping.sh`

```bash
#!/bin/bash
#SBATCH --job-name=mapping_bwamem
#SBATCH -p normal
#SBATCH -c 8
#SBATCH --array=0-88%4
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=4G

RAW_DATA="/scratch/MOryzae/DATA/Raw_Data"
REF_GENOME="/scratch/MOryzae/Reference/MOryzae_70-15.fna"
BAM_OUTPUT="/scratch/MOryzae/Results/Mapping"
mkdir -p ${BAM_OUTPUT}

# Get sample name from array index
SAMPLES=($(ls ${RAW_DATA}/*_R1.fastq.gz | xargs basename -a | sed 's/_R1.fastq.gz//'))
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

module load bwa/0.7.17
module load samtools/1.15

# Perform alignment
bwa mem -t 8 \
        -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPL:ILLUMINA" \
        ${REF_GENOME} \
        ${RAW_DATA}/${SAMPLE}_R1.fastq.gz \
        ${RAW_DATA}/${SAMPLE}_R2.fastq.gz | \
samtools sort -@ 8 -o ${BAM_OUTPUT}/${SAMPLE}.sorted.bam

# Index the BAM file
samtools index ${BAM_OUTPUT}/${SAMPLE}.sorted.bam

echo "Mapping completed for: ${SAMPLE}"
```

**Submission:**
```bash
sbatch Mapping.sh
```

---

### 4. Variant Calling

#### SNP and Indel Discovery (GATK HaplotypeCaller)

**SLURM Script:** `VariantCalling.sh`

```bash
#!/bin/bash
#SBATCH --job-name=variant_calling
#SBATCH -p normal
#SBATCH -c 4
#SBATCH --array=0-88%8
#SBATCH --time=08:00:00
#SBATCH --mem-per-cpu=4G

BAM_DIR="/scratch/MOryzae/Results/Mapping"
REF_GENOME="/scratch/MOryzae/Reference/MOryzae_70-15.fna"
VCF_OUTPUT="/scratch/MOryzae/Results/SNP/individual_gvcfs"
mkdir -p ${VCF_OUTPUT}

SAMPLES=($(ls ${BAM_DIR}/*.sorted.bam | xargs basename -a | sed 's/.sorted.bam//'))
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

module load gatk/4.3.0
module load samtools

# Call variants per-sample (GVCF format)
gatk HaplotypeCaller \
     -R ${REF_GENOME} \
     -I ${BAM_DIR}/${SAMPLE}.sorted.bam \
     -O ${VCF_OUTPUT}/${SAMPLE}.g.vcf.gz \
     -ERC GVCF

echo "GVCF created for: ${SAMPLE}"
```

#### Joint Genotyping (GATK GenomicsDBImport & GenotypeGVCFs)

```bash
#!/bin/bash
#SBATCH --job-name=joint_genotyping
#SBATCH -p normal
#SBATCH -c 8
#SBATCH --time=12:00:00
#SBATCH --mem=64G

GVCF_DIR="/scratch/MOryzae/Results/SNP/individual_gvcfs"
REF_GENOME="/scratch/MOryzae/Reference/MOryzae_70-15.fna"
DB_DIR="/scratch/MOryzae/Results/SNP/genomicsdb"
VCF_JOINT="/scratch/MOryzae/Results/SNP/joint_calls.vcf.gz"

module load gatk/4.3.0

# Create list of GVCFs
ls ${GVCF_DIR}/*.g.vcf.gz > gvcf_list.txt

# Import GVCFs into GenomicsDB
gatk GenomicsDBImport \
     --sample-name-map gvcf_list.txt \
     --genomics-db-workspace ${DB_DIR} \
     -L all_contigs.list \
     --reader-threads 8

# Perform joint genotyping
gatk GenotypeGVCFs \
     -R ${REF_GENOME} \
     -V gendb://${DB_DIR} \
     -O ${VCF_JOINT} \
     --allow-old-rms-mapping-quality-annotation-data

echo "Joint genotyping completed. Output: ${VCF_JOINT}"
```

---

### 5. Variant Filtering and Quality Control

#### SNP Quality Filtering

**Rationale:** Remove low-confidence variants that inflate Type I error in downstream population genetics inferences.

**SLURM Script:** `FilterVariants.sh`

```bash
#!/bin/bash
#SBATCH --job-name=filter_vcf
#SBATCH -p normal
#SBATCH -c 4
#SBATCH --time=04:00:00
#SBATCH --mem=16G

VCF_INPUT="/scratch/MOryzae/Results/SNP/joint_calls.vcf.gz"
VCF_FILTERED="/scratch/MOryzae/Results/SNP/vcf_filtered/filtered_snps.vcf"
QC_OUTPUT="/scratch/MOryzae/Results/SNP/qc_stats"
mkdir -p ${QC_OUTPUT} $(dirname ${VCF_FILTERED})

module load gatk/4.3.0
module load vcftools

# Filtering thresholds (standard for population genetics):
# - QUAL ≥ 30
# - DP per-sample ≥ 5, ≤ 100
# - MQ ≥ 40
# - Missing data ≤ 20%

gatk SelectVariants \
     -V ${VCF_INPUT} \
     --select-type-to-include SNP \
     -O ${QC_OUTPUT}/snps_only.vcf.gz

gatk VariantFiltration \
     -V ${QC_OUTPUT}/snps_only.vcf.gz \
     --filter-expression "QUAL < 30" --filter-name "LowQual" \
     --filter-expression "QD < 2.0" --filter-name "LowQD" \
     --filter-expression "MQ < 40" --filter-name "LowMQ" \
     --filter-expression "FS > 60" --filter-name "HighFS" \
     -O ${QC_OUTPUT}/snps_filtered.vcf.gz

# Apply filters (remove FAIL variants)
gatk SelectVariants \
     -V ${QC_OUTPUT}/snps_filtered.vcf.gz \
     --exclude-filtered \
     -O ${VCF_FILTERED}.gz

# Filter by missing data (max 20% missing per SNP)
vcftools --gzvcf ${VCF_FILTERED}.gz \
         --max-missing 0.8 \
         --recode \
         --out ${VCF_FILTERED%.*}

echo "Filtering completed. Variants retained: $(zcat ${VCF_FILTERED}.gz | grep -v '^#' | wc -l)"
```

---

### 6. Population Genetics Analyses

The following analyses were performed on the filtered SNP dataset to investigate population structure, genetic diversity, and signatures of divergence.

#### 6.1 Genetic Diversity Metrics

**Objective:** Quantify within-population genetic variation.

**R Script:** `genetic_diversity.R`

```r
#!/usr/bin/env Rscript

# Load required libraries
library(pegas)          # Population genetics analysis
library(ggplot2)        # Data visualization
library(dplyr)          # Data manipulation

# Set working directory
setwd("/scratch/MOryzae/Results/Population_Structure")

# Read filtered VCF
vcf_file <- "/scratch/MOryzae/Results/SNP/vcf_filtered/filtered_snps.vcf"
vcf_data <- read.vcf(vcf_file)

# Extract genotype matrix
genotypes <- extract.gt(vcf_data)
genotypes_numeric <- as.numeric(genotypes)

# Calculate diversity indices per population
compute_diversity_metrics <- function(genotypes, pop_assignments) {
  
  diversity_results <- list()
  
  for (pop in unique(pop_assignments)) {
    pop_genotypes <- genotypes[, pop_assignments == pop]
    
    # Nucleotide diversity (π)
    pi <- nucleotide.div(pop_genotypes, model = "raw")
    
    # Expected heterozygosity (He)
    he <- expected.het(pop_genotypes)
    
    # Observed heterozygosity (Ho)
    ho <- observed.het(pop_genotypes)
    
    # Number of alleles
    n_alleles <- colSums(pop_genotypes > 0)
    
    diversity_results[[pop]] <- data.frame(
      Population = pop,
      Nucleotide_Diversity = mean(pi, na.rm = TRUE),
      Expected_Heterozygosity = mean(he, na.rm = TRUE),
      Observed_Heterozygosity = mean(ho, na.rm = TRUE),
      Mean_Alleles = mean(n_alleles, na.rm = TRUE)
    )
  }
  
  do.call(rbind, diversity_results)
}

# Assign populations based on prior knowledge or clustering results
pop_assignments <- factor(rep(c("rice", "wheat", "grass"), each = 30))

diversity_table <- compute_diversity_metrics(genotypes, pop_assignments)
write.csv(diversity_table, "genetic_diversity_by_population.csv", row.names = FALSE)

print(diversity_table)
```

#### 6.2 Population Structure: Model-Based Clustering (LEA/SNMF)

**Objective:** Infer the most likely number of genetic clusters (K) and estimate admixture proportions.

**R Script:** `population_structure_lea.R`

```r
#!/usr/bin/env Rscript

library(LEA)
library(vcfR)
library(ggplot2)
library(tidyr)

# Set working directory
lea_dir <- "/scratch/MOryzae/Results/Population_Structure/LEA"
dir.create(lea_dir, recursive = TRUE, showWarnings = FALSE)

# Read filtered VCF
vcf_file <- "/scratch/MOryzae/Results/SNP/vcf_filtered/filtered_snps.vcf"
vcf_data <- read.vcfR(vcf_file)

# Extract sample IDs
sample_ids <- colnames(vcf_data@gt)[-1]

# Convert VCF to GENO format (required by LEA)
cat("Converting VCF to GENO format...\n")
output <- vcf2geno(vcf_file, output.file = file.path(lea_dir, "filtered_snps"))

# Define GENO file path
geno_file <- file.path(lea_dir, "filtered_snps.geno")

# Run SNMF for K = 1 to 10 (test range for K)
cat("Running SNMF for K = 1 to 10...\n")
project <- snmf(geno_file,
                K = 1:10,
                entropy = TRUE,
                repetitions = 10,
                project = "new",
                seed = 12345)

# Extract cross-entropy values
entropy_values <- cross.entropy(project, K = 1:10)
entropy_df <- data.frame(K = 1:10, Cross_Entropy = entropy_values)

# Identify optimal K
optimal_k <- which.min(entropy_values)
cat(sprintf("Optimal K = %d (lowest cross-entropy = %.4f)\n", optimal_k, min(entropy_values)))

# Plot entropy criterion
png(file.path(lea_dir, "entropy_criterion.png"), width = 1600, height = 800)
plot(project, col = "steelblue", pch = 19, cex = 2,
     xlab = "Number of Clusters (K)",
     ylab = "Cross-Entropy",
     main = "LEA/SNMF Cross-Entropy Criterion")
abline(v = optimal_k, col = "red", lty = 2, lwd = 2)
text(optimal_k, max(entropy_values) * 0.95, sprintf("K = %d", optimal_k), 
     pos = 4, col = "red", font = 2)
dev.off()

# Extract ancestry matrix (Q-matrix) for optimal K
best_run <- which.min(cross.entropy(project, K = optimal_k))
ancestry_matrix <- Q(project, K = optimal_k, run = best_run)

# Add sample IDs to ancestry matrix
ancestry_with_ids <- data.frame(
  Sample_ID = sample_ids,
  ancestry_matrix
)
colnames(ancestry_with_ids)[2:(optimal_k+1)] <- paste0("Cluster_", 1:optimal_k)

write.csv(ancestry_with_ids, 
          file.path(lea_dir, sprintf("ancestry_matrix_K%d.csv", optimal_k)),
          row.names = FALSE)

cat(sprintf("Ancestry matrix saved for K = %d\n", optimal_k))

# Visualize ancestry matrix (barplot)
# Sort individuals by dominant cluster
dominant_cluster <- apply(ancestry_matrix, 1, which.max)
sort_order <- order(dominant_cluster, decreasing = TRUE)
ancestry_sorted <- ancestry_matrix[sort_order, ]
sample_ids_sorted <- sample_ids[sort_order]

# Create barplot
png(file.path(lea_dir, sprintf("ancestry_barplot_K%d.png", optimal_k)), 
    width = 2000, height = 900)

barplot(t(ancestry_sorted),
        col = rainbow(optimal_k),
        border = NA,
        space = 0,
        xlab = "Individual",
        ylab = "Ancestry Proportion",
        main = sprintf("Population Structure (K = %d)", optimal_k),
        cex.axis = 0.8)

# Add sample labels
axis(1, at = 1:length(sample_ids_sorted), labels = sample_ids_sorted,
     las = 2, cex.axis = 0.6)

dev.off()

cat(sprintf("Ancestry barplot saved: %s\n", 
            file.path(lea_dir, sprintf("ancestry_barplot_K%d.png", optimal_k))))

# Create cluster assignment table
cluster_assignments <- data.frame(
  Sample_ID = sample_ids,
  Dominant_Cluster = dominant_cluster
)

write.csv(cluster_assignments,
          file.path(lea_dir, "cluster_assignments.csv"),
          row.names = FALSE)
```

#### 6.3 Phylogenetic Reconstruction

**Objective:** Estimate evolutionary relationships and assess divergence times between populations.

**Workflow:**

```bash
# Convert VCF to PHYLIP format
vcftools --vcf filtered_snps.vcf \
         --out filtered_snps \
         --phylip

# Infer phylogeny using RAxML-NG
raxml-ng --parse \
         --msa filtered_snps.phy \
         --model GTR+G+ASC_FELSENSTEIN \
         --prefix alignment_check

raxml-ng --search \
         --msa filtered_snps.phy \
         --model GTR+G+ASC_FELSENSTEIN \
         --prefix MOryzae_tree \
         --threads 8 \
         --seed 12345

# Perform 1000 bootstrap replicates for support
raxml-ng --bootstrap \
         --msa filtered_snps.phy \
         --model GTR+G+ASC_FELSENSTEIN \
         --prefix MOryzae_bootstrap \
         --threads 8 \
         --num-replicates 1000

# Map bootstrap values onto best tree
raxml-ng --support \
         --tree MOryzae_tree.raxml.bestTree \
         --bs-trees MOryzae_bootstrap.raxml.bootstraps \
         --prefix MOryzae_final
```

**Visualization (R):**

```r
library(ape)
library(ggtree)

tree <- read.tree("MOryzae_final.raxml.supportTree")

# Root tree at midpoint
tree <- midpoint.root(tree)

# Plot with bootstrap support values
ggtree(tree, aes(color = bootstrap), size = 1.2) +
  geom_tiplab(size = 4) +
  geom_nodepoint(aes(size = bootstrap), alpha = 0.6) +
  scale_color_gradient(low = "lightgray", high = "red") +
  theme_tree() +
  ggtitle("Phylogenetic Tree of M. oryzae Isolates with Bootstrap Support")

ggsave("phylogeny_with_bootstrap.png", width = 14, height = 10, dpi = 300)
```

#### 6.4 Linkage Disequilibrium Assessment

```r
library(LDcorSV)
library(ggplot2)

# Calculate pairwise LD (r²) for SNPs within genomic windows
genotypes_numeric <- apply(genotypes, 2, as.numeric)

ld_matrix <- ld(genotypes_numeric, depth = 5000)

# Summarize LD decay
ld_summary <- data.frame(
  Distance_kb = seq_len(nrow(ld_matrix)) * 1,
  Mean_r_squared = rowMeans(ld_matrix, na.rm = TRUE)
)

ggplot(ld_summary, aes(x = Distance_kb, y = Mean_r_squared)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", color = "blue", fill = "lightblue") +
  labs(x = "Genomic Distance (kb)", y = "Mean r²", 
       title = "Linkage Disequilibrium Decay") +
  theme_minimal()

ggsave("ld_decay.png", width = 10, height = 6, dpi = 300)
```

---

## Results

Comprehensive results from all analyses are organized in the [Results](/Results) directory:

```
Results/
├── Mapping/                      # BAM files and alignment statistics
├── SNP/
│   ├── vcf_filtered/            # Filtered SNP VCF file
│   └── qc_stats/                # Variant filtering statistics
├── Population_Structure/
│   ├── LEA/                     # LEA/SNMF ancestry matrices, barplots
│   └── Diversity_Metrics/       # Genetic diversity summary tables
├── Phylogenetics/               # Phylogenetic trees with bootstrap support
└── Visualizations/              # Publication-ready figures
```

### Key Findings Summary

[Placeholder: Results to be populated with empirical findings]

---

## Project Organization and Documentation

### Repository Structure

```
Mentoring_Project/
├── README.md                          # This file
├── Topic/
│   ├── Sujet.pdf                     # Project specification (French)
│   └── Data_project_pyri.xlsx        # Metadata and sample information
├── Jupyter_books/                     # Jupyter notebooks documenting analyses
│   ├── 01_QualityControl.ipynb
│   ├── 02_ReadMapping.ipynb
│   ├── 03_VariantCalling.ipynb
│   ├── 04_PopulationStructure.ipynb
│   ├── 05_PhylogeneticInference.ipynb
│   └── 06_DataVisualization.ipynb
├── Wrappers/                          # Standalone R scripts for analysis
│   ├── genetic_diversity.R
│   ├── population_structure_lea.R
│   └── phylogenetic_visualization.R
├── SCRIPTS/                           # HPC submission scripts
│   ├── FastQC.sh
│   ├── Mapping.sh
│   ├── VariantCalling.sh
│   ├── FilterVariants.sh
│   └── [additional scripts]
├── Reference/                         # Reference genome and associated files
│   └── MOryzae_70-15.fna
└── Results/                           # Analysis outputs (see above)
```

### Mind Map and Analysis Roadmap

A comprehensive mind map summarizing the analytical strategy, methodological considerations, and key decision points is available at:  
[Interactive Mind Map](https://mm.tt/app/map/3489792035?t=sMvnsWUbpm)

The mind map encompasses:
- **Thematic overview** of population genetics concepts
- **Literature synthesis** on *M. oryzae* biology and genetic diversity
- **Bioinformatic strategy** flow chart
- **Analytical workflow** with software tools and parameters
- **Results interpretation** framework

---

## Methods Summary

### Software and Versions

| Tool | Version | Purpose | Reference |
|------|---------|---------|-----------|
| FastQC | 0.11.9+ | Read quality assessment | Andrews (2010) |
| MultiQC | 1.9+ | Aggregated QC reporting | Ewels et al. (2016) |
| BWA-MEM | 0.7.17+ | Read mapping | Li & Durbin (2009) |
| SAMtools | 1.15+ | BAM manipulation | Li et al. (2009) |
| GATK | 4.3.0+ | Variant calling & filtering | McKenna et al. (2010) |
| VCFtools | 0.1.16+ | VCF processing | Danecek et al. (2011) |
| LEA | 3.0+ | Population structure inference | Frichot & François (2015) |
| RAxML-NG | 1.0+ | Phylogenetic inference | Kozlov et al. (2019) |
| R | 4.0+ | Statistical analysis | R Core Team (2021) |

### Statistical Methods

**Population Structure Inference:** Latent Factor Mixed Models (LEA) via the SNMF algorithm, which applies non-negative matrix factorization to uncover hidden population structure. Model selection was performed using the cross-entropy criterion to identify the optimal number of clusters.

**Phylogenetic Reconstruction:** Maximum Likelihood inference under the General Time-Reversible (GTR) substitution model with gamma-distributed rate heterogeneity (Γ). Bootstrap resampling (1,000 replicates) was used to assess topological support.

**Genetic Diversity:** Standard population genetics metrics including nucleotide diversity (π), expected heterozygosity, and observed heterozygosity were computed using the *pegas* R package.

---

## Reproducibility

### System Requirements

- **Operating System:** Linux (Ubuntu 20.04+ recommended)
- **HPC Scheduler:** SLURM
- **Programming Languages:** Bash, R (≥4.0)
- **Memory:** 64 GB minimum for joint genotyping
- **Storage:** ~1 TB for raw reads + intermediate files

### Data Availability

Raw sequencing data is archived under IRD project governance. Processed data (filtered VCF, ancestry matrices) are available upon request to the project mentors or through the corresponding GitHub repository.

### Code Availability

All scripts, Jupyter notebooks, and analysis code are freely available in this repository and conform to open science principles.

### Computational Reproducibility

To reproduce analyses on your own system:

1. Install required software (see Software and Versions table)
2. Download the [Singularity/Docker container](link-if-available) containing all dependencies, OR manually install tools listed in the table
3. Obtain filtered VCF file from the [Results](/Results) directory
4. Execute scripts from the [Wrappers](/Wrappers) directory in numerical order
5. Outputs will be written to a `Results/` subdirectory

**Example:**
```bash
Rscript Wrappers/genetic_diversity.R
Rscript Wrappers/population_structure_lea.R
```

---

## Learning Outcomes

Participants in this project developed proficiency in:

- **Bioinformatics Workflow Design:** Planning and executing multi-stage computational pipelines for genomic data
- **High-Performance Computing:** SLURM job submission, resource optimization, and batch processing on HPC clusters
- **Quality Assurance:** Implementing rigorous quality control checks throughout the analysis pipeline
- **Population Genomics:** Applying model-based and phylogenetic approaches to study population structure
- **Scientific Communication:** Documenting analyses with sufficient detail for independent reproduction
- **Collaboration:** Working within an international research network on complex biological problems

---

## References

Andrews, S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data. *Babraham Institute*, UK.

Danecek, P., Auton, A., Abecasis, G., et al. (2011). The variant call format and VCFtools. *Bioinformatics*, 27(15), 2156–2158. https://doi.org/10.1093/bioinformatics/btr330

Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: Summarize analysis results for multiple tools and samples in a single report. *Bioinformatics*, 32(19), 3047–3048. https://doi.org/10.1093/bioinformatics/btw354

Frichot, E., & François, O. (2015). LEA: An R package for landscape and ecological association studies. *Methods in Ecology and Evolution*, 6(8), 925–929. https://doi.org/10.1111/2041-210X.12382

Kozlov, A. M., Darriba, D., Flouri, T., Morel, B., & Stamatakis, A. (2019). RAxML-NG: A fast, scalable, and user-friendly tool for phylogenetic analysis. *Bioinformatics*, 35(21), 4453–4455. https://doi.org/10.1093/bioinformatics/btz305

Li, H., & Durbin, R. (2009). Fast and accurate short read alignment with Burrows-Wheeler Transform. *Bioinformatics*, 25(14), 1754–1760. https://doi.org/10.1093/bioinformatics/btp324

Li, H., Handsaker, B., Wysoker, A., et al. (2009). The Sequence Alignment/Map format and SAMtools. *Bioinformatics*, 25(16), 2078–2079. https://doi.org/10.1093/bioinformatics/btp352

McKenna, A., Hanna, M., Banks, E., et al. (2010). The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data. *Genome Biology*, 11, R47. https://doi.org/10.1186/gb-2010-11-9-r47

R Core Team. (2021). *R: A Language and Environment for Statistical Computing*. R Foundation for Statistical Computing, Vienna, Austria.

---

## Citation

If you use this project's code, methods, or results, please cite:

```bibtex
@misc{name_magnaporthe_2024,
  author = {NAME, Pakyendou Estel and AGNIMONHAN, Attolou Raoul},
  title = {Population Genetics and Cryptic Species Assessment in \textit{Magnaporthe oryzae}},
  year = {2024},
  howpublished = {GitHub Repository},
  url = {https://github.com/estelp/Mentoring_Project},
  note = {Training project, International Certificate in Bioinformatics and Genomics (CIBiG)}
}
```

---

## License

This project is licensed under the [specify license: MIT, CC-BY-4.0, GPL-3.0]. The analysis workflows and scripts are freely available for educational and research purposes.

---

## Acknowledgments

We gratefully acknowledge:

- **Aurore COMTE** and **Sébastien RAVEL** for expert mentorship and methodological guidance
- The **WAVE Initiative** for computational infrastructure and project support
- **AfricaRice** and **INERA** for institutional support
- The **IRD Montpellier Cluster** for high-performance computing resources

---

## Frequently Asked Questions (FAQ)

### Q1: Can I reuse these scripts on my own *M. oryzae* samples?

**A:** Yes, with modifications. The scripts are parameterized and can be adapted to different datasets by changing path variables (e.g., `RAW_DATA_PATH`, sample array limits, reference genome). However, ensure your sample metadata, sequencing platform, and organism are compatible with the chosen reference genome.

### Q2: What sequencing depth is required for robust SNP calling?

**A:** For population-level inference, we recommend a minimum of 10× genome coverage per sample. This project used higher coverage to maximize SNP discovery. If lower coverage data is available, adjust the depth filtering parameters in `FilterVariants.sh`.

### Q3: How long does the pipeline take to complete?

**A:** For 89 samples on the Montpellier cluster with 8-core allocations:
- Quality Control: ~2 hours
- Read Mapping: ~12–15 hours (8 samples/batch)
- Variant Calling: ~24 hours (joint genotyping step is rate-limiting)
- Filtering & Analysis: ~4–6 hours

Total: ~2–3 days wall-clock time.

### Q4: Where can I find the input data?

**A:** Raw sequencing data is stored on the IRD cluster (restricted access). Filtered VCF files and results are provided in the [Results](/Results) directory.

### Q5: What if my variants fail to pass filtering thresholds?

**A:** Adjust filtering criteria in `FilterVariants.sh` if:
- Data are lower coverage: increase `--max-missing` threshold
- Specific population has high indel rates: modify `FS` (Fisher Strand) threshold
- Quality scores are systematically lower: reduce `QUAL` minimum threshold

Always verify changes using variant effect predictors (SnpEff, VEP) to ensure functional variants are retained.

---

## Contributing

We welcome contributions, bug reports, and suggestions for pipeline improvements. Please open an issue or submit a pull request on the [GitHub repository](https://github.com/estelp/Mentoring_Project).

---

**Last Updated:** November 2024  
**Corresponding Contact:** Pakyendou Estel NAME (pakyendou.name@gmail.com)
