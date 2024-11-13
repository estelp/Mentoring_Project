#!/bin/bash

############# SLURM Configuration ##############

### Set the job name
#SBATCH --job-name=mapping_pipeline

### Set the partition to use
#SBATCH -p normal

### Set the number of CPUs to use
#SBATCH -c 16

### Specify the node on which the job should run
#SBATCH --nodelist=node20  # Specifies that the job should run on node20

#################################################

########### Path Variables ##################

# Define paths for data and output directories
REF_PATH="/scratch/MOryzae/REF/MOryzae_genomic.fna"
TRIM_DATA_PATH="/scratch/MOryzae/DATA/Trimming"
OUTPUT_PATH="/scratch/MOryzae/MAPPING"
SAM_PATH="/scratch/MOryzae/MAPPING/sam_files"
BAM_PATH="/scratch/MOryzae/MAPPING/bam_raw"
STATS_PATH="/scratch/MOryzae/MAPPING/bam_stats"
FILTERED_PATH="/scratch/MOryzae/MAPPING/bam_filtered"
SORTED_PATH="/scratch/MOryzae/MAPPING/bam_mapped_sort"
STAT_FILE="${STATS_PATH}/all_stat.csv"

# Load bwa-mem2 and samtools modules
module load bwamem2/2.2.1
module load samtools/1.18

# List of sequences
sequences=("12-1-205" "BN0019" "Br80" "CHRF" "GN0001" "IN0114" "JP0091" "Pg1213-22" "PY86-1" 
           "TN0002" "WBKY11" "87-120" "BN0119" "CD0065" "CHW" "GRF52" "IN0115" "LpKY97" 
           "PgKY" "SSFL02" "TN0050" "WHTQ" "AG0004" "BN0123" "CD0142" "CM0028" "GY0040" 
           "IN0116" "ML0060" "PGPA" "SSFL14-3" "TN0057" "Y34" "Arcadia" "BN0202" "CD0156" 
           "EI9411" "HO" "INA168" "ML0062" "PH42" "SV9610" "TN0065" "Z2-1" "B2" "BN0252" 
           "CD0156" "EI9604" "IA1" "IR00102" "ML33" "PL2-1" "SV9623" "TN0090" "B51" 
           "BR0032" "CH0043" "FH" "IB33" "IR0013" "NG0012" "PL3-1" "T25" "TR0025" "B71" 
           "BR0032" "CH0072" "FR1067" "IB49" "IR0015" "NG0054" "PY0925" "TF05-1" "US0041" 
           "Bd8401" "Br130" "CH0452" "FR1069" "IC17" "IR0083" "NP0058" "PY36-1" "TG0004" 
           "US0064" "BdBar" "Br48" "CH0461" "G17" "IE1K" "IR0084" "P131" "PY5033" "TG0032" 
           "US0071" "BdMeh" "Br58" "CH0533" "G22" "IN0017" "IR0088" "P28" "PY5033" 
           "TH0016" "US0071" "BF0072" "Br62" "CH1103" "GFSI1-7-2" "IN0054" "IR0095" 
           "P29" "PY6017" "TH0016" "VT0027" "Bm88324" "Br7" "CH1164" "GG11" "IN0059" 
           "IT0010" "P3" "PY6045" "TN0001" "VT0030")

#################################################

########### Create output directories ##################

# Create the output directories if they do not exist
mkdir -p "$SAM_PATH" "$BAM_PATH" "$STATS_PATH" "$FILTERED_PATH" "$SORTED_PATH"

########### Step 1: Mapping, SAM to BAM conversion, Statistics, Filtering, Sorting ##################

# Create the CSV statistics file
echo "sample,mapped,primary_mapped,properly_paired,unmapped" > "$STAT_FILE"

# Loop over each sequence to perform all steps
for sequence in "${sequences[@]}"; do
    echo -e "######################\nProcessing for ${sequence}..."
    
    # Define file paths for input and output
    R1="${TRIM_DATA_PATH}/${sequence}_R1_paired.fastq.gz"
    R2="${TRIM_DATA_PATH}/${sequence}_R2_paired.fastq.gz"
    SAM_FILE="${SAM_PATH}/${sequence}.sam"
    BAM_FILE="${BAM_PATH}/${sequence}.bam"
    FLAGSTAT_FILE="${STATS_PATH}/${sequence}.flagstat"
    FILTERED_BAM="${FILTERED_PATH}/${sequence}.mappedpaired.bam"
    SORTED_BAM="${SORTED_PATH}/${sequence}.mappedpaired.sorted.bam"
    SORTED_BAM_INDEX="${SORTED_PATH}/${sequence}.mappedpaired.sorted.bam.bai"
    
    # Check if input files exist
    if [[ -f "$R1" && -f "$R2" ]]; then
        # Step 1: Mapping with bwa-mem2
        bwa-mem2 mem -t 16 "$REF_PATH" "$R1" "$R2" -o "$SAM_FILE"
        if [[ $? -ne 0 ]]; then
            echo "Mapping error for ${sequence}. Skipping."
            continue
        fi
        echo "Mapping completed for ${sequence}"
    else
        echo "Error: Missing files for ${sequence}. Check ${R1} and ${R2}."
        continue
    fi

    # Step 2: Convert SAM to BAM
    if [[ -f "$SAM_FILE" ]]; then
        samtools view -b -o "$BAM_FILE" "$SAM_FILE"
        if [[ $? -ne 0 ]]; then
            echo "Error during SAM to BAM conversion for ${sequence}. Skipping."
            continue
        fi
        echo "SAM to BAM conversion successful for ${sequence}"
    else
        echo "Error: Missing SAM file for ${sequence}."
        continue
    fi

    # Step 3: Generate statistics using flagstat
    if [[ -f "$BAM_FILE" ]]; then
        samtools flagstat -@ 16 "$BAM_FILE" > "$FLAGSTAT_FILE"
        if [[ $? -ne 0 ]]; then
            echo "Error extracting statistics for ${sequence}. Skipping."
            continue
        fi
        echo "Statistics generated for ${sequence}: ${FLAGSTAT_FILE}"

        # Extract statistics data and append it to the CSV file
        mapped=$(grep "mapped (" "$FLAGSTAT_FILE" | awk '{print $1}')
        primary_mapped=$(grep "primary paired (" "$FLAGSTAT_FILE" | awk '{print $1}')
        properly_paired=$(grep "properly paired (" "$FLAGSTAT_FILE" | awk '{print $1}')
        unmapped=$(grep "unmapped (" "$FLAGSTAT_FILE" | awk '{print $1}')
        
        echo "${sequence},${mapped},${primary_mapped},${properly_paired},${unmapped}" >> "$STAT_FILE"
    else
        echo "Error: Missing BAM file for ${sequence}."
        continue
    fi

    # Step 4: Filter BAM files
    if [[ -f "$BAM_FILE" ]]; then
        samtools view -bh -@ 16 -f 0x02 -o "$FILTERED_BAM" "$BAM_FILE"
        if [[ $? -ne 0 ]]; then
            echo "Error during BAM filtering for ${sequence}. Skipping."
            continue
        fi
        echo "Filtered BAM created for ${sequence}: ${FILTERED_BAM}"
    else
        echo "Error: Missing BAM file for ${sequence}."
        continue
    fi

    # Step 5: Sort the filtered BAM files
    if [[ -f "$FILTERED_BAM" ]]; then
        samtools sort -@ 16 "$FILTERED_BAM" -o "$SORTED_BAM"
        if [[ $? -ne 0 ]]; then
            echo "Error during BAM sorting for ${sequence}. Skipping."
            continue
        fi
        echo "Sorted BAM created for ${sequence}: ${SORTED_BAM}"
    else
        echo "Error: Missing filtered BAM file for ${sequence}."
        continue
    fi

    # Step 6: Index the sorted BAM file
    if [[ -f "$SORTED_BAM" ]]; then
        samtools index "$SORTED_BAM"
        if [[ $? -ne 0 ]]; then
            echo "Error during BAM indexing for ${sequence}. Skipping."
            continue
        fi
