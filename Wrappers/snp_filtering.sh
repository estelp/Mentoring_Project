#!/bin/bash

#SBATCH --job-name=snp_filtering
#SBATCH -p normal
#SBATCH -c 8
#SBATCH --nodelist=node20

# Définir les répertoires
INPUT_DIR="/scratch/MOryzae/SNP/vcf_files"
OUTPUT_DIR="/scratch/MOryzae/SNP/vcf_filtered"

# Charger les modules nécessaires
module load bcftools/1.18
module load vcftools/0.1.16
module load htslib/1.19

# Paramètres de filtrage
MAF=0.1
MISS=0.9
QUAL=19000

# Fichiers
VCF_IN="${INPUT_DIR}/snp_correct.vcf.gz"
VCF_OUT="${OUTPUT_DIR}/filtered_snps.vcf.gz"
VCF_STATS="${OUTPUT_DIR}/filtered_snps_stats.txt"

# Vérifier que le répertoire de sortie existe
if [ ! -d "$OUTPUT_DIR" ]; then
    mkdir -p "$OUTPUT_DIR"
    echo "Répertoire de sortie créé : $OUTPUT_DIR"
fi

# Vérifier que le fichier d'entrée existe
if [ ! -f "$VCF_IN" ]; then
    echo "Erreur : Le fichier d'entrée $VCF_IN n'existe pas."
    exit 1
fi

# Exécuter vcftools pour filtrer les SNPs
vcftools --gzvcf $VCF_IN \
    --remove-indels \
    --maf $MAF \
    --max-missing $MISS \
    --minQ $QUAL \
    --recode \
    --stdout | bgzip -c > $VCF_OUT

# Vérifier que le fichier de sortie a été créé
if [ ! -f "$VCF_OUT" ]; then
    echo "Erreur : Le fichier $VCF_OUT n'a pas été généré."
    exit 1
fi
echo "Filtrage terminé. Fichier filtré disponible à : $VCF_OUT"

# Indexation avec bcftools
bcftools index $VCF_OUT
if [ $? -ne 0 ]; then
    echo "Erreur lors de l'indexation du fichier filtré."
    exit 1
fi
echo "Indexation terminée."

# Calcul des statistiques avec bcftools
bcftools stats $VCF_OUT > $VCF_STATS
if [ $? -ne 0 ]; then
    echo "Erreur lors du calcul des statistiques."
    exit 1
fi
echo "Statistiques disponibles dans : $VCF_STATS"

