#!/bin/bash

############# SLURM Configuration ##############

### Define Job name
#SBATCH --job-name=genome_download_index

### Define partition to use
#SBATCH -p normal

### Define number of CPUs to use
#SBATCH -c 8

### Specify the node to run on
#SBATCH --nodelist=node20  # Spécifie que le job doit être exécuté sur node20

#################################################

########### Execution Command ###################

# Créer le répertoire REF si nécessaire
REF_DIR="/scratch/MOryzae/REF"
mkdir -p "$REF_DIR"

# Définir le chemin du fichier de référence
GENOME_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/495/GCF_000002495.2_MG8/GCF_000002495.2_MG8_genomic.fna.gz"
GENOME_FILE="$REF_DIR/MOryzae_genomic.fna.gz"

# Télécharger le génome de référence
echo "Téléchargement du génome de référence..."
wget -O "$GENOME_FILE" "$GENOME_URL"

# Vérifier si le téléchargement a réussi
if [[ $? -ne 0 ]]; then
    echo "Erreur : Le téléchargement du génome a échoué."
    exit 1
fi

# Décompresser le fichier
echo "Décompression du fichier..."
gunzip "$GENOME_FILE"

# Modifier le nom pour une utilisation plus facile
mv "${GENOME_FILE%.gz}" "$REF_DIR/MOryzae_genomic.fna"

# Charger le module bwa-mem2 pour faire l'indexation
module load bwamem2/2.2.1

# Indexer le génome de référence
echo "Indexation du génome de référence..."
bwa-mem2 index "$REF_DIR/MOryzae_genomic.fna"

# Vérifier si l'indexation a réussi
if [[ $? -eq 0 ]]; then
    echo "Indexation du génome de référence réussie."
else
    echo "Erreur : L'indexation du génome a échoué."
fi
