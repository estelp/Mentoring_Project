#!/bin/bash

############# SLURM Configuration ##############

### Set the job name
#SBATCH --job-name=mapping_pipeline

### Set the partition to use
#SBATCH -p normal

### Set the number of CPUs to use
#SBATCH -c 8

### Specify the node on which the job should run
#SBATCH --nodelist=node20  # Specifies that the job should run on node20

#################################################

########### Path Variables ##################

# Define paths for data and output directories
REF_PATH="/scratch/MOryzae/REF/MOryzae_genomic.fna"
TRIM_DATA_PATH="/scratch/MOryzae/DATA/Trimming"
OUTPUT_PATH="/scratch/MOryzae/MAPPING"
SAM_PATH="${OUTPUT_PATH}/sam_files"
BAM_PATH="${OUTPUT_PATH}/bam_raw"
STATS_PATH="${OUTPUT_PATH}/bam_stats"
FILTERED_PATH="${OUTPUT_PATH}/bam_filtered"
SORTED_PATH="${OUTPUT_PATH}/bam_mapped_sort"
STAT_FILE="${STATS_PATH}/all_stat.csv"

# Load bwa-mem2 and samtools modules
module load bwamem2/2.2.1
module load samtools/1.18

# List of sequences
sequences=("AG0004" "BN0123" "CH0461" "G22" "IE1K" "IR0015" "ML33" "PH42" "TN0057"
           "Arcadia" "BN0202" "CH0533" "GFSI1-7-2" "IN0017" "IR0083" "NG0012" "PL2-1" "TN0065"
           "B2" "BN0252" "CH1103" "GG11" "IN0054" "IR0084" "NG0054" "SSFL02" "TN0090"
           "B71" "Br7" "CH1164" "GN0001" "IN0059" "IR0088" "NP0058" "SSFL14-3" "TR0025"
           "Bd8401" "Br80" "CHRF" "GY0040" "IN0114" "IR0095" "P28" "T25" "US0041"
           "BdBar" "CD0065" "CHW" "HO" "IN0115" "IT0010" "P29" "TG0004" "US0064"
           "BF0072" "CD0142" "CM0028" "IA1" "IN0116" "JP0091" "P3" "TG0032" "VT0027"
           "Bm88324" "CH0043" "FR1067" "IB33" "INA168" "LpKY-97-1" "Pg1213-22" "TN0001" "VT0030"
           "BN0019" "CH0072" "FR1069" "IB49" "IR00102" "ML0060" "PgKY4OV2-1" "TN0002" "Z2-1"
           "BN0119" "CH0452" "G17" "IC17" "IR0013" "ML0062" "PgPA18C-02" "TN0050")

#################################################

########### Create output directories ##################

# Create the output directories if they do not exist
mkdir -p "$SAM_PATH" "$BAM_PATH" "$STATS_PATH" "$FILTERED_PATH" "$SORTED_PATH"

########### Step 1: Mapping, SAM to BAM conversion, Statistics, Filtering, Sorting ##################

# Create the CSV statistics file and add headers
echo "sample,total_reads,total_reads_perc,mapped_reads,mapped_reads_perc,properly_paired,properly_paired_perc,singletons,singletons_perc,unmapped_reads,unmapped_reads_perc,read1,read1_perc,read2,read2_perc" > "$STAT_FILE"

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
    
    # Step 1: Mapping with bwa-mem2
    bwa-mem2 mem -t 8 "$REF_PATH" "$R1" "$R2" -o "$SAM_FILE"
    echo "Mapping completed for ${sequence}"
    
    # Step 2: Convert SAM to BAM
    samtools view -b -o "$BAM_FILE" "$SAM_FILE"
    echo "SAM to BAM conversion successful for ${sequence}"
    
    # Step 3: Generate statistics using flagstat
    samtools flagstat -@ 8 "$BAM_FILE" > "$FLAGSTAT_FILE"
    echo "Statistics generated for ${sequence}: ${FLAGSTAT_FILE}"

    # Step 4: Filter BAM files
    samtools view -bh -@ 8 -f 0x02 -o "$FILTERED_BAM" "$BAM_FILE"
    echo "Filtered BAM created for ${sequence}: ${FILTERED_BAM}"
    
    # Step 5: Sort the filtered BAM files
    samtools sort -@ 8 "$FILTERED_BAM" -o "$SORTED_BAM"
    echo "Sorted BAM created for ${sequence}: ${SORTED_BAM}"
    
    # Step 6: Index the sorted BAM file
    samtools index "$SORTED_BAM"
    echo "Indexing completed for ${sequence}"
    
done

# Process each flagstat file and write results to CSV with percentages
for file in "$STATS_PATH"/*flagstat; do
    if [[ -f $file ]]; then
        sample_name=$(basename "$file" | cut -d. -f1)
        new_line="$sample_name,"
        
        # Extract values from the flagstat file
        total_reads=$(grep "in total" "$file" | awk '{print $1}')
        mapped_reads=$(grep "mapped (" "$file" | awk '{print $1}')
        mapped_reads_perc=$(grep "mapped (" "$file" | awk -F '[()%]' '{print $2}')
        
        properly_paired=$(grep "properly paired (" "$file" | awk '{print $1}')
        properly_paired_perc=$(grep "properly paired (" "$file" | awk -F '[()%]' '{print $2}')
        
        singletons=$(grep "singletons (" "$file" | awk '{print $1}')
        singletons_perc=$(grep "singletons (" "$file" | awk -F '[()%]' '{print $2}')
        
        unmapped_reads=$((total_reads - mapped_reads))
        unmapped_reads_perc=$(awk "BEGIN {print (100 - $mapped_reads_perc)}")
        
        read1=$(grep "read1" "$file" | awk '{print $1}')
        read1_perc=$(awk "BEGIN {print ($read1/$total_reads)*100}")
        
        read2=$(grep "read2" "$file" | awk '{print $1}')
        read2_perc=$(awk "BEGIN {print ($read2/$total_reads)*100}")
        
        # Append statistics with percentages to the CSV file
        new_line+="${total_reads},100,${mapped_reads},${mapped_reads_perc},${properly_paired},${properly_paired_perc},${singletons},${singletons_perc},${unmapped_reads},${unmapped_reads_perc},${read1},${read1_perc},${read2},${read2_perc}"
        echo "$new_line" >> "$STAT_FILE"
    fi
done
