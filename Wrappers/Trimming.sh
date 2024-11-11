#!/bin/bash

############# SLURM Configuration ##############

### Define Job name
#SBATCH --job-name=trimmomatic

### Define partition to use
#SBATCH -p normal

### Define number of CPUs to use
#SBATCH -c 8

### Specify the node to run on
#SBATCH --nodelist=node20  # Spécifie que le job doit être exécuté sur node20

### Define array job (ajustez la plage selon le nombre de fichiers)
#SBATCH --array=0-88%4  # 89 fichiers (index de 0 à 88)

#################################################

########### Execution Command ##################

# Charger le module Trimmomatic
module load Trimmomatic/0.39

# Définir les chemins vers les répertoires de travail
RAW_DATA_PATH="/scratch/MOryzae/DATA/Raw_Data"
TRIMMOMATIC_OUTPUT_PATH="/scratch/MOryzae/DATA/Trimming"

# Liste des fichiers sans suffixe
files=("AG0004" "BN0123" "CH0461" "G22" "IE1K" "IR0015" "ML33" "PH42" "TN0057"
       "Arcadia" "BN0202" "CH0533" "GFSI1-7-2" "IN0017" "IR0083" "NG0012" "PL2-1" "TN0065"
       "B2" "BN0252" "CH1103" "GG11" "IN0054" "IR0084" "NG0054" "SSFL02" "TN0090"
       "B71" "Br7" "CH1164" "GN0001" "IN0059" "IR0088" "NP0058" "SSFL14-3" "TR0025"
       "Bd8401" "Br80" "CHRF" "GY0040" "IN0114" "IR0095" "P28" "T25" "US0041"
       "BdBar" "CD0065" "CHW" "HO" "IN0115" "IT0010" "P29" "TG0004" "US0064"
       "BF0072" "CD0142" "CM0028" "IA1" "IN0116" "JP0091" "P3" "TG0032" "VT0027"
       "Bm88324" "CH0043" "FR1067" "IB33" "INA168" "LpKY-97-1" "Pg1213-22" "TN0001" "VT0030"
       "BN0019" "CH0072" "FR1069" "IB49" "IR00102" "ML0060" "PgKY4OV2-1" "TN0002" "Z2-1"
       "BN0119" "CH0452" "G17" "IC17" "IR0013" "ML0062" "PgPA18C-02" "TN0050")

# Récupérer l'index du job SLURM
index=$SLURM_ARRAY_TASK_ID

# Définir les noms de fichiers en fonction de l'index du tableau
file1="${files[$index]}_R1.fastq.gz"
file2="${files[$index]}_R2.fastq.gz"

# Créer le répertoire de sortie s'il n'existe pas
mkdir -p "$TRIMMOMATIC_OUTPUT_PATH"


# Utiliser Trimmomatic pour le trimming des reads avec un score Phred >= 30
echo "Running Trimmomatic on $file1 and $file2..."

trimmomatic PE -phred33 "$RAW_DATA_PATH/$file1" "$RAW_DATA_PATH/$file2" \
    "$TRIMMOMATIC_OUTPUT_PATH/${files[$index]}_R1_paired.fastq.gz" \
    "$TRIMMOMATIC_OUTPUT_PATH/${files[$index]}_R1_unpaired.fastq.gz" \
    "$TRIMMOMATIC_OUTPUT_PATH/${files[$index]}_R2_paired.fastq.gz" \
    "$TRIMMOMATIC_OUTPUT_PATH/${files[$index]}_R2_unpaired.fastq.gz" \
    SLIDINGWINDOW:4:30 \
    LEADING:3 TRAILING:3 \
    MINLEN:36

# Vérifier si Trimmomatic a réussi
if [[ $? -eq 0 ]]; then
    echo "Trimmomatic completed successfully for $file1 and $file2."
else
    echo "Erreur : Trimmomatic a rencontré un problème avec $file1 et $file2."
    exit 1
fi
