#!/bin/bash

############# SLURM Configuration ##############

### Define Job name
#SBATCH --job-name=multiqc

### Define partition to use
#SBATCH -p normal

### Define number of CPUs to use
#SBATCH -c 8

### Specify the node to run on
#SBATCH --nodelist=node20  # Spécifie que le job doit être exécuté sur node20

#################################################

########### Execution Command ##################

# Définir les chemins vers les répertoires de travail
FASTQC_PATH="/scratch/MOryzae/QC/FastQC_Trim"
MULTIQC_OUTPUT_PATH="/scratch/MOryzae/QC/MultiQC_Trim"

# Charger le module MultiQC
module load multiqc/1.9

# Créer le répertoire de sortie s'il n'existe pas
mkdir -p "$MULTIQC_OUTPUT_PATH"

# Vérifier si des rapports FastQC existent dans le répertoire d'entrée
if [ ! -d "$FASTQC_PATH" ] || [ -z "$(ls -A $FASTQC_PATH)" ]; then
    echo "Error: No FastQC reports found in $FASTQC_PATH."
    exit 1  # Quitte avec un code d'erreur si aucun rapport n'est trouvé
fi

# Lancer MultiQC sur les rapports de FastQC
echo "Running MultiQC on FastQC reports in $FASTQC_PATH..."
multiqc "$FASTQC_PATH" -o "$MULTIQC_OUTPUT_PATH"

# Vérifier si MultiQC a réussi
if [[ $? -eq 0 ]]; then
    echo "MultiQC completed successfully."
else
    echo "Error: MultiQC encountered an issue."
    exit 1  # Quitter avec un code d'erreur si MultiQC échoue
fi
