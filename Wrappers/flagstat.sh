#!/bin/bash

############# SLURM Configuration ##############

### Set the job name
#SBATCH --job-name=concatenate_flagstat

### Set the partition to use
#SBATCH -p normal

### Set the number of CPUs to use
#SBATCH -c 8

### Specify the node on which the job should run
#SBATCH --nodelist=node20  # Specifies that the job should run on node20

#################################################

# Path variables
flagstat_dir="/scratch/MOryzae/MAPPING/bam_stats"  # Directory containing flagstat files
stat_file="/scratch/MOryzae/MAPPING/bam_stats/all_stat.csv"  # Output file

# List of prefix
prefix=("AG0004" "BN0123" "CH0461" "G22" "IE1K" "IR0015" "ML33" "PH42" "TN0057"
           "Arcadia" "BN0202" "CH0533" "GFSI1-7-2" "IN0017" "IR0083" "NG0012" "PL2-1" "TN0065"
           "B2" "BN0252" "CH1103" "GG11" "IN0054" "IR0084" "NG0054" "SSFL02" "TN0090"
           "B71" "Br7" "CH1164" "GN0001" "IN0059" "IR0088" "NP0058" "SSFL14-3" "TR0025"
           "Bd8401" "Br80" "CHRF" "GY0040" "IN0114" "IR0095" "P28" "T25" "US0041"
           "BdBar" "CD0065" "CHW" "HO" "IN0115" "IT0010" "P29" "TG0004" "US0064"
           "BF0072" "CD0142" "CM0028" "IA1" "IN0116" "JP0091" "P3" "TG0032" "VT0027"
           "Bm88324" "CH0043" "FR1067" "IB33" "INA168" "LpKY-97-1" "Pg1213-22" "TN0001" "VT0030"
           "BN0019" "CH0072" "FR1069" "IB49" "IR00102" "ML0060" "PgKY4OV2-1" "TN0002" "Z2-1"
           "BN0119" "CH0452" "G17" "IC17" "IR0013" "ML0062" "PgPA18C-02" "TN0050")

# Create the output file and write the headers
echo "Sequence,total,primary,secondary,supplementary,duplicates,primary_duplicates,mapped,primary_mapped,paired_in_sequencing,read1,read2,properly_paired,singletons,unmapped" > "$stat_file"

# Loop through each sequence
for file in "${prefix[@]}"; do
    # Define the corresponding flagstat file for the sequence
    flagstat_file="${flagstat_dir}/${file}.flagstat"

    # Check if the file exists
    if [[ -f "$flagstat_file" ]]; then
        # Initialize variables to store values from the flagstat file
        total=0
        primary=0
        secondary=0
        supplementary=0
        duplicates=0
        primary_duplicates=0
        mapped=0
        primary_mapped=0
        paired_in_sequencing=0
        read1=0
        read2=0
        properly_paired=0
        singletons=0
        unmapped=0

        # Read the flagstat file line by line
        while IFS= read -r line; do
            # Extract values based on the specific keywords
            if [[ $line =~ ^([0-9]+)\ \+\ [0-9]+\ in\ total ]]; then
                total="${BASH_REMATCH[1]}"
            elif [[ $line =~ ^([0-9]+)\ \+\ [0-9]+\ primary ]]; then
                primary="${BASH_REMATCH[1]}"
            elif [[ $line =~ ^([0-9]+)\ \+\ [0-9]+\ secondary ]]; then
                secondary="${BASH_REMATCH[1]}"
            elif [[ $line =~ ^([0-9]+)\ \+\ [0-9]+\ supplementary ]]; then
                supplementary="${BASH_REMATCH[1]}"
            elif [[ $line =~ ^([0-9]+)\ \+\ [0-9]+\ duplicates ]]; then
                duplicates="${BASH_REMATCH[1]}"
            elif [[ $line =~ ^([0-9]+)\ \+\ [0-9]+\ primary\ duplicates ]]; then
                primary_duplicates="${BASH_REMATCH[1]}"
            elif [[ $line =~ ^([0-9]+)\ \+\ [0-9]+\ mapped ]]; then
                mapped="${BASH_REMATCH[1]}"
            elif [[ $line =~ ^([0-9]+)\ \+\ [0-9]+\ primary\ mapped ]]; then
                primary_mapped="${BASH_REMATCH[1]}"
            elif [[ $line =~ ^([0-9]+)\ \+\ [0-9]+\ paired\ in\ sequencing ]]; then
                paired_in_sequencing="${BASH_REMATCH[1]}"
            elif [[ $line =~ ^([0-9]+)\ \+\ [0-9]+\ read1 ]]; then
                read1="${BASH_REMATCH[1]}"
            elif [[ $line =~ ^([0-9]+)\ \+\ [0-9]+\ read2 ]]; then
                read2="${BASH_REMATCH[1]}"
            elif [[ $line =~ ^([0-9]+)\ \+\ [0-9]+\ properly\ paired ]]; then
                properly_paired="${BASH_REMATCH[1]}"
            elif [[ $line =~ ^([0-9]+)\ \+\ [0-9]+\ singletons ]]; then
                singletons="${BASH_REMATCH[1]}"
            fi
        done < "$flagstat_file"

        # Calculate unmapped reads
        unmapped=$((total - mapped))

        # Write the results to the CSV file
        echo "$file,$total,$primary,$secondary,$supplementary,$duplicates,$primary_duplicates,$mapped,$primary_mapped,$paired_in_sequencing,$read1,$read2,$properly_paired,$singletons,$unmapped" >> "$stat_file"
    else
        echo "Warning: File for sequence $seq not found. Skipping."
    fi
done

echo "Data extraction complete. Results saved in $stat_file."
