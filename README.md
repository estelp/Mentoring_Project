# Mentoring_Project

# Bioinformatics Tutored Project

## Overview

This README provides comprehensive information about a practical training course in Bioinformatics and Genomics, aligned with the International Certificate in Bioinformatics and Genomics (CIBiG). The course takes place in Abidjan, Côte d'Ivoire, as part of the Central and West African Virus Epidemiology (WAVE) initiative.

## Project Details

- **Project Duration:** October 21, 2024 – November 25, 2024
- **Course Dates:** September 9, 2024 – October 4, 2024
- **Location:** Abidjan, Côte d'Ivoire

## Partners

This project is a collaborative effort involving:

- **AGNIMONHAN Attolou Raoul** (AfricaRice, Ivory Coast)
- **NAME Pakyendou Estel** (INERA-WAVE, Burkina Faso)

## Tutors

The project will be guided by experienced professionals:

- **Aurore COMTE** (IRD)
- **Sébastien RAVEL** (CIRAD)

## Objectives

The primary aim of this tutored project is to provide hands-on training in bioinformatics techniques and applications, focusing on genomic data analysis and interpretation, particularly in the context of viral epidemiology. Participants will engage in practical exercises, collaborative projects, and discussions led by the tutors.

## Expected Outcomes

Participants will gain:
- Proficiency in bioinformatics tools and software.
- Skills in genomic data analysis.
- Understanding of the role of bioinformatics.
- Collaboration experience with professionals in the field.

## Contact Information

For further inquiries or information about the project, please contact the tutors or project partners:

- **Aurore COMTE:** [aurore.comte@ird.fr]
- **Sébastien RAVEL:** [sebastien.ravel@cirad.fr]
- **AGNIMONHAN Attolou Raoul:** [R.Agnimonhan@cgiar.org]
- **NAME Pakyendou Estel:** [pakyendou.name@gmail.com]


## Project Topic

[Access the directory](/Topic/Sujet.pdf)

<u>Specific lineages of <em>Magnaporthe oryzae</em></u>

Studying the genetic structures of pathogen populations, in relation to life-history traits such as mode of reproduction, host range or resistance to treatments, is essential for understanding the emergence and spread of infectious diseases. Among plant pathogens, the ascomycete fungus <em>Magnaporthe oryzae</em>, responsible for blast disease in many cultivated and wild grass species, is a model of interest. Although this pathogen is mainly studied for its devastating effects on rice (<em>Oryza sativa</em>), it also infects other cereal crops, such as wheat, barley and millet, as well as grasses such as ryegrass and St Augustine's grass. Previous research has shown that <em>M. oryzae</em> is subdivided into several host-specific lineages, with genetic divergence probably linked to host changes.
This study aims to further understand the genetic structure of several Magnaporthe isolates ([see table](/Topic/Data_project_pyri.xlsx)), from different host species, to determine whether they form host-specific lineages and to assess the existence of cryptic species within <em>M. oryzae</em>.
In a few words:
We are trying to understand the population structure of <em>M. oryzae</em>.
● What is the relationship between host and population structure?
● Are there cryptic species that stand out from the rest of the <em>M. oryzae</em> population, or is <em>M. oryzae</em> made up of just one species, independent of its host?


## Mind Map

Acces(https://mm.tt/app/map/3489792035?t=sMvnsWUbpm)

This mind map summarizes in one image all the reflections made around our theme. It enabled us to better understand the subject we were given, then define the different approaches to analysis and rendering. It includes:
- Thematic
- Bibliographic synthesis
- Databases 
- Methodology
- Bioinformatics strategy
- Results
- Restitution & Reproducibility


## BIOINFORMTIC STRATEGY

The entire bioinformatics strategy was carried out in several stages. 
The jupyters for all these stages are available in the [Jupyter_books](/Jupyter_books) directory.

### 1. Data acquisition

All work or bioinformatics strategy was carried out on the Montpellier cluster (bioinfo-master1.ird.fr).
So, the first step before executing the scripts was 

* to connect to the cluster on your terminal via your login

```bash
ssh login@bioinfo-master1.ird.fr
```

* Connect to a node on a specific partition according to your analyses
x = partition
y = number of cpus

```bash
srun -p x -c y --pty bash -i
```

* Create the “MOryzae” working directory in the “scratch” directory and move around in it

```bash
mkdir -p /scratch/MOryzae
cd /scratch/MOryzae
```

* Display absolute path

```bash
pwd
```

* Display directory contents

```bash
ls -lh
```

* Create a “DATA” subdirectory for sequencing data

```bash
mkdir DATA
```

* Create two sub-directories “Raw_Data” for paired-end sequences and “Contigs” for assembled sequences 

```bash
mkdir -p DATA/Raw_Data 
mkdir -p DATA/Contigs
```

* Copy the sequencing data from the NAS into the “DATA/Raw_Data” directory

```bash
scp -r san:/projects/medium/CIBiG_MOryzae/fastq/*_R* /scratch/MOryzae/DATA/Raw_Data
```
* List Raw_Data content and view one sequence

```bash
ls -lh /scratch/MOryzae/DATA/Raw_Data/
head /scratch/MOryzae/DATA/Raw_Data/x.fastq.gz  
### x = sequence name
```

### 2. Quality Control

#### 2.1. Analysis

This step is crucial in order to judge the quality of the reads and to know what steps to take before proceeding with the actual analyses.

Create a “QC” directory in the working directory in which quality control results will be stored and Move to the “QC” directory
mkdir -p /scratch/MOryzae/QC

```bash
mkdir -p /scratch/MOryzae/QC
cd /scratch/MOryzae/QC
```
* **<u>FASTQC</u>**

Create the “FastQC” directory in the QC directory for outputs 

```bash
mkdir FastQC
```

Create the SCRIPTS directory in the working directory
(this directory will be used to store all scripts to be run in sbatch mode on the cluster)

```bash
mkdir -p /scratch/MOryzae/SCRIPTS
```
Move to the created directory

```bash
cd /scratch/MOryzae/SCRIPTS
```

Open nano text editor

```bash
nano FastQC.sh
```
save the following sbatch script

```bash
#!/bin/bash

############# SLURM Configuration ##############

### Define Job name
#SBATCH --job-name=qc_fastq

### Define partition to use
#SBATCH -p normal

### Define number of CPUs to use
#SBATCH -c 2

### Define array job (ajustez la plage selon le nombre de fichiers)
#SBATCH --array=0-88%4  # 89 paires de fichiers (0 à 88)

### Specify the node to run on
#SBATCH --nodelist=node20  # Spécifie que le job doit être exécuté sur node20

#################################################

########### Execution Command ##################

# Define working and output directory absolute path
RAW_DATA_PATH="/scratch/MOryzae/DATA/Raw_Data"
QC_OUTPUT_PATH="/scratch/MOryzae/QC/FastQC"

# Load  "FastQC" module available
module load FastQC/0.11.9

# Create output directory

mkdir -p "$QC_OUTPUT_PATH"


# List of files without suffixe
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

# Get file index 
index=$SLURM_ARRAY_TASK_ID

# Get file name from index
file1="${files[$index]}_R1.fastq.gz"
file2="${files[$index]}_R2.fastq.gz"

# Check if file existence in directory before run the command
if [[ -f "$RAW_DATA_PATH/$file1" && -f "$RAW_DATA_PATH/$file2" ]]; then
    echo "Running FastQC on $file1 and $file2..."
    fastqc "$RAW_DATA_PATH/$file1" "$RAW_DATA_PATH/$file2" -o "$QC_OUTPUT_PATH"

    # Check if command succes
    if [[ $? -eq 0 ]]; then
        echo "FastQC completed successfully for $file1 and $file2."
    else
        echo "Error: FastQC encountered an issue with $file1 and $file2."
    fi
else
    echo "Error: One or both files do not exist: $file1, $file2."
    exit 1  # Quitter avec un code d'erreur
fi


```

Run the script 
[Access FastqQC.sh](/Wrappers/FastQC.sh)

```bash
sbash FastQC.sh
```

At the end of the task, check the contents

```bash
ls -lh /scratch/MOryzae/QC/FastQC/
```

Use generated html files to check read quality


* **<u>MULTIQC</u>**

This step consolidates the html reports from FastQC into a single, easy-to-interpret report. To do this

Create the “MultiQC” directory in the QC directory for outputs 


```bash
mkdir -p /scratch/MOryzae/QC/MultiQC
```

Open the nano text editor to edit a sbatch script

```bash
nano MultiQC.sh
```

save the following sbatch script

```bash
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
FASTQC_PATH="/scratch/MOryzae/QC/FastQC"
MULTIQC_OUTPUT_PATH="/scratch/MOryzae/QC/MultiQC"

# Charger le module MultiQC
module load multiqc/1.9

# Créer le répertoire de sortie

mkdir -p "$MULTIQC_OUTPUT_PATH"


# Lancer MultiQC sur les rapports de FastQC
echo "Running MultiQC on FastQC reports in $FASTQC_PATH..."
multiqc "$FASTQC_PATH" -o "$MULTIQC_OUTPUT_PATH"

# Vérifier si MultiQC a réussi
if [[ $? -eq 0 ]]; then
    echo "MultiQC completed successfully."
else
    echo "Error: MultiQC encountered an issue."
fi


```

Run the script
[Access MultiQC.sh](/Wrappers/MultiQC.sh)

```bash
sbash MultiQC.sh
```

At the end of the task, check the contents

```bash
ls -lh /scratch/MOryzae/QC/MultiQC/
```

Use generated html files to check read quality

For the next step,
Move the entire contents of the QC directory to the NAS

```bash
scp -r /scratch/MOryzae/QC san:/projects/medium/CIBiG_MOryzae/
```

Retrieve this QC directory from the NAS on your local machine to analyze the results

```bash
scp -r login@bioinfo-san.ird.fr:/projects/medium/CIBiG_MOryzae/QC /path/to/working dorectory/on your laptop/
```

### 2.2. Trimming

After analyzing the quality control results, we noticed that some sequences were not of good quality.
As a result, trimming of the reads was necessary. 
The Trimmomatic tool will be used, setting the Phred score at 30 to retain only the best quality reads for further analysis.

* **<u>TRIMMING</u>**

Create a “Trimming” subdirectory in the DATA directory

```bash
mkdir -p /scratch/MOryzae/DATA/Trimming
```

Make sure you're in the /SCRIPTS directory

```bash
cd /scratch/MOryzae/SCRIPTS/
```

Open nano text editor

```bash
nano Trimming.sh
```

save the following sbatch script

```bash
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

```

Run the script
[Access Trimming.sh](/Wrappers/Trimming.sh)

```bash
sbash Trimming.sh
```

At the end of the task, check the contents

```bash
ls -lh /scratch/MOryzae/DATA/Trimming/
```

restart quality control on trimmed sequences


* **<u>FASTQC_TRIM</u>**

Create the “FastQC_Trim” directory in the QC directory for outputs

```bash
mkdir -p /scratch/MOryzae/DATA/FastQC_Trim
```

Create the SCRIPTS directory in the working directory
(this directory will be used to store all scripts to be run in sbatch mode on the cluster)


```bash
cd /scratch/MOryzae/SCRIPTS
```

Open nano text editor

```bash
nano FastQC_Trim.sh
```

save the following sbatch script

```bash
#!/bin/bash

############# SLURM Configuration ##############

### Define Job name
#SBATCH --job-name=qc_fastq

### Define partition to use
#SBATCH -p normal

### Define number of CPUs to use
#SBATCH -c 2

### Define array job (ajustez la plage selon le nombre de fichiers)
#SBATCH --array=0-88%4  # 89 paires de fichiers (0 à 88)

### Specify the node to run on
#SBATCH --nodelist=node20  # Spécifie que le job doit être exécuté sur node20

#################################################

########### Execution Command ##################

# Define working and output directory absolute path
TRIM_DATA_PATH="/scratch/MOryzae/DATA/Trimming"
QC_OUTPUT_PATH="/scratch/MOryzae/QC/FastQC_Trim"

# Load "FastQC" module available
module load FastQC/0.11.9

# Create output directory if it doesn't exist
mkdir -p "$QC_OUTPUT_PATH"

# List of files without suffix
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

# Get file index
index=$SLURM_ARRAY_TASK_ID

# Get file name from index
file1="${files[$index]}_R1_paired.fastq.gz"
file2="${files[$index]}_R2_paired.fastq.gz"

# Check if files exist before running FastQC
if [[ -f "$TRIM_DATA_PATH/$file1" && -f "$TRIM_DATA_PATH/$file2" ]]; then
    echo "Running FastQC on $file1 and $file2..."

    # Run FastQC
    fastqc "$TRIM_DATA_PATH/$file1" "$TRIM_DATA_PATH/$file2" -o "$QC_OUTPUT_PATH"

    # Check if FastQC ran successfully
    if [[ $? -eq 0 ]]; then
        echo "FastQC completed successfully for $file1 and $file2."
    else
        echo "Error: FastQC encountered an issue with $file1 and $file2."
        exit 1  # Exit with error code if FastQC fails
    fi
else
    echo "Error: One or both files do not exist: $file1, $file2."
    exit 1  # Exit with error code if files are missing
fi

```

Run the script
[Access FastqQC_Trim.sh](/Wrappers/FastQC_Trim.sh)

```bash
sbash FastQC_Trim.sh
```

At the end of the task, check the contents

```bash
ls -lh /scratch/MOryzae/QC/FastQC_Trim/
```

Use generated html files to check read quality


* **<u>MULTIQC_TRIM</u>**

This step consolidates the html reports from FastQC into a single, easy-to-interpret report. To do this

Create the “MultiQC_Trim” directory in the QC directory for outputs

```bash
mkdir -p /scratch/MOryzae/QC/MultiQC_Trim
```

Open the nano text editor to edit a sbatch script

```bash
nano MultiQC_Trim.sh
```

save the following sbatch script


```bash
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

```

Run the script
[Access MultiQC_Trim.sh](/Wrappers/MultiQC_Trim.sh)

```bash
sbash MultiQC_Trim.sh
```

At the end of the task, check the contents

```bash
ls -lh /scratch/MOryzae/QC/MultiQC_Trim/
```

Use generated html files to check read quality

For the next step,
Move the entire contents of the QC directory to the NAS

```bash
scp -r /scratch/MOryzae/QC san:/projects/medium/CIBiG_MOryzae/
```

Retrieve this QC directory from the NAS on your local machine to analyze the results

```bash
scp -r login@bioinfo-san.ird.fr:/projects/medium/CIBiG_MOryzae/QC /path/to/working dorectory/on your laptop/
```


### 3. Mapping

### 3.1. Reference genome download & indexing

In order to go a little faster, we've written a sbatch script for this step. 
The script takes into account reference genome download, decompression and indexing.

Move to the SCRIPTS directory

```bash
cd /scratch/MOryzae/SCRIPTS
```

Open nano text editor

```bash
nano Refseq.sh
```

save the following sbatch script

```bash
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

```

Run the script
[Access Refseq.sh](/Wrappers/Refseq.sh)

```bash
sbash Refseq.sh
```

At the end of the task, check the contents

```bash
ls -lh /scratch/MOryzae/REF/
```

### 3.2. Mapping Run

This script includes not only the mapping command, but also commands for:
. converting sam files to bam
. mapping statistics
. grouping statistics into a csv file
. filtering of mapped reads
. output of bam files containing mapped reads
. and indexing of these bam files

Move to the SCRIPTS directory

```bash
cd /scratch/MOryzae/SCRIPTS
```

Open nano text editor

```bash
nano mapping_pipeline.sh
```

save the following sbatch script

```bash
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

```

Run the script [Access mapping_pipeline.sh](/Wrappers/mapping_pipeline.sh)

```bash
sbash mapping_pipeline.sh
```

Group flagstat files into a single csv file

Open nano text editor

```bash
nano flagstat.sh
```

save the following sbatch script

```bash
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

```

Run the script
[Access flagstat.sh](/Wrappers/flagstat.sh)

```bash
sbash flagstat.sh
```


## Delete the sam_files, bam_raw and bam_filtered subdirectories in the MAPPING directory to free up more space.

```bash
rm -rf /scratch/MOryzae/MAPPING/sam_files
rm -rf /scratch/MOryzae/MAPPING/bam_raw
rm -rf /scratch/MOryzae/MAPPING/bam_filtered
```

At the end of the task, check the contents

```bash
ls -lh /scratch/MOryzae/MAPPING
```

Use generated csv file to interpret this mapping step

For the next step,
Move the entire contents of the MAPPING directory to the NAS

```bash
scp -r /scratch/MOryzae/MAPPING san:/projects/medium/CIBiG_MOryzae/
```

Retrieve this MAPPING directory from the NAS on your local machine to analyze the results

```bash
scp -r login@bioinfo-san.ird.fr:/projects/medium/CIBiG_MOryzae/MAPPING /path/to/working dorectory/on your laptop/
```

### 4. SNP calling

SNP calling (detection of SNPs, or Single Nucleotide Polymorphisms) is a fundamental analysis in genomics and bioinformatics, aimed at identifying genetic variations of a single nucleotide in the genome of an organism.
In our context, this analysis was carried out for:
... Studying genetic diversity: SNPs are common genetic markers and provide information on variation between individuals or populations. This is essential for understanding genetic diversity, evolution and the specific adaptations of populations.
... Evolution and phylogeny: SNP analysis can be used to trace evolutionary relationships between species and to establish phylogenetic trees, helping to understand evolutionary history and the links between populations or species.

### 4.1. Reference genome indexing

Move to the reference directory

```bash
cd /scratch/MOryzae/REF
```

Load samtools module 

```bash
module load samtools/1.18
```

Reference indexing using samtools faidx

```bash
samtools faidx /scratch/MOryzae/REF/MOryzae_genomic.fna
```

Check directory contents and generated files

```bash
ls -lh
```

### 4.2. Run SNP Calling


Move to the SCRIPTS directory

```bash
cd /scratch/MOryzae/SCRIPTS
```

Open nano text editor

```bash
nano snp_calling_pipeline.sh
```

save the following sbatch script

```bash
#!/bin/bash

############ SLURM Configuration ##############

### Define Job name
#SBATCH --job-name=snp2_calling

### Define partition to use
#SBATCH -p normal

### Define number of CPUs to use
#SBATCH -c 16

### Specify the node to run on
#SBATCH --nodelist=node20  # Run the job on node20

#################################################

########### Execution Commands ###################

# Variables
SORTED_PATH="/scratch/MOryzae/MAPPING/bam_mapped_sort"
VCF_OUTPUT="/scratch/MOryzae/SNP/vcf_files/all_samples.vcf.gz"
SNP_OUTPUT="/scratch/MOryzae/SNP/vcf_files/snp.vcf.gz"
REF_GENOME="/scratch/MOryzae/REF/MOryzae_genomic.fna"
SNP_STATS_DIR="/scratch/MOryzae/SNP/stats"

# Load necessary modules
module load samtools/1.18
module load bcftools/1.18

# Create necessary directories
mkdir -p /scratch/MOryzae/SNP/vcf_files "$SNP_STATS_DIR"

# Step 1: Generate VCF file, compress, and index
echo -e "######################\nGenerating compressed VCF file"
bcftools mpileup -Ou --threads 16 -f "$REF_GENOME" "$SORTED_PATH"/*.mappedpaired.sorted.bam | \
  bcftools call -mv -Oz -o "$VCF_OUTPUT" || {
    echo "Error: VCF generation and compression failed" >&2
    exit 1
}

# Step 2: Index the compressed VCF file
echo -e "######################\nIndexing compressed VCF file"
bcftools index "$VCF_OUTPUT" || {
    echo "Error: VCF file indexing failed" >&2
    exit 1
}

# Step 3: Filter to retain only SNPs
echo -e "######################\nFiltering to retain only SNPs"
bcftools view -v snps -Oz -o "$SNP_OUTPUT" "$VCF_OUTPUT" || {
    echo "Error: SNP filtering failed" >&2
    exit 1
}

# Step 4: Index the SNP VCF file
echo -e "######################\nIndexing SNP VCF file"
bcftools index "$SNP_OUTPUT" || {
    echo "Error: SNP VCF file indexing failed" >&2
    exit 1
}

# Step 5: Generate SNP statistics
echo -e "######################\nGenerating SNP statistics"
bcftools stats "$SNP_OUTPUT" > "$SNP_STATS_DIR/all_samples_SNP_statistics.txt" || {
    echo "Error: Failed to generate SNP statistics" >&2
    exit 1
}

echo "Compressed VCF file and SNP-specific VCF file generated and indexed successfully."


```

Run the script
[Access snp_calling_pipeline.sh](/Wrappers/snp_calling_pipeline.sh)

```bash
sbash snp_calling_pipeline.sh
```

At the end of the task, check the contents

```bash
ls -lh /scratch/MOryzae/SNP
```

Move to the SNP directory

```bash
cd /scratch/MOryzae/SNP/vcf_files
```

Load bcftools module

```bash
module load bcftools/1.18
```

Use bcftools to list the names of current samples:

```bash
bcftools query -l snp.vcf.gz> snp_samples.txt
```

Modify samples.txt or create a new file, e.g. new_samples.txt, 
mapping the old names to the new simplified names.

```bash
awk -F'/' '{print $0 "\t" $NF}' snp_samples.txt | sed 's/.mappedpaired.sorted.bam//g' > new_snp_samples.txt
```

Recover isolate names only in a new text file

```bash
awk -F'/' '{print $0 "\t" $NF}' new_snp_samples.txt | cut -f3 > new2_snp_samples.txt
```

Use bcftools reheader to apply the changes:

```bash
bcftools reheader -s new2_snp_samples.txt -o snp_correct.vcf.gz snp.vcf.gz
```

To check that the new names have been applied correctly:

```bash
bcftools query -l snp_correct.vcf.gz
```

Use generated files to interpret this step and plannig next analysis

For the next step,
Move the entire contents of the SNP directory to the NAS

```bash
scp -r /scratch/MOryzae/SNP san:/projects/medium/CIBiG_MOryzae/
```

Retrieve this SNP directory from the NAS on your local machine to analyze the results

```bash
scp -r login@bioinfo-san.ird.fr:/projects/medium/CIBiG_MOryzae/SNP /path/to/working dorectory/on your laptop/
```

### 4.3. PCA & DAPC

#### Evaluate level of missing data (by sample, by positions)

Create a PLINK directory in the working directory

```bash
mkdir -p /scratch/MOryzae/PLINK
```

Move to the created directory

```bash
cd PLINK/
```

Load plink module first

```bash
module load plink/1.9
```

Start evaluation of missing data

```bash
plink -vcf /scratch/MOryzae/SNP/vcf_files/snp_correct.vcf.gz --allow-extra-chr --missing --out  ./plink/dataset
```

#### Generate PCA using genotyping information contained in VCF

Plink alllows to create a PCA (principal components analysis) of samples, so that we can easily evaluate genetic distance between samples.

This will generate a matrix of coordinates in the different component. By default, it provides the first 20 principal components of the variance-standardized relationship matrix. We will focus only the first 3 axes for subsequent visualization (--pca 3)


Run following commands

```bash
plink -vcf /scratch/MOryzae/SNP/vcf_files/snp_correct.vcf.gz --allow-extra-chr --cluster --matrix --pca 3 --mind --out ./plink/dataset
```

#### Convert "eigenvec" generate file to csv file

Move to the SCRIPTS directory

```bash
cd /scratch/MOryzae/SCRIPTS
```

Open nano text editor

```bash
nano eigenvec_to_csv.sh
```

save the following sbatch script

```bash
#!/bin/bash

############# SLURM Configuration ##############

### Define Job name
#SBATCH --job-name=eigenvec_to_csv

### Define partition to use
#SBATCH -p normal

### Define number of CPUs to use
#SBATCH -c 8

### Specify the node to run on
#SBATCH --nodelist=node20  # Run the job on node20

#################################################

########### Execution Command ###################

# Define directories
INPUT_DIR="/scratch/MOryzae/PLINK"
OUTPUT_DIR="/scratch/MOryzae/PLINK"

# List of directories to process
DIRECTORIES=("plink")

# Loop through each directory
for DIR in "${DIRECTORIES[@]}"; do
    PCA_FILE="$INPUT_DIR/$DIR/dataset.eigenvec"
    OUTPUT_CSV="$OUTPUT_DIR/$DIR/dataset.csv"

    # Check if PCA results exist
    if [ -f "$PCA_FILE" ]; then
        echo "Processing PCA results for $DIR..."

        # Convert eigenvec file to CSV format
        awk 'NR==1{print "FID,IID,PC1,PC2,PC3"} NR>1{print $1","$2","$3","$4","$5}' "$PCA_FILE" > "$OUTPUT_CSV"

        # Check if conversion was successful
        if [ $? -eq 0 ]; then
            echo "Conversion successful: $OUTPUT_CSV created."
        else
            echo "Error: Failed to convert $PCA_FILE to CSV."
            exit 1
        fi
    else
        echo "Error: PCA results file not found in $DIR."
        exit 1
    fi
done

echo "PCA analysis completed. Results are saved in $OUTPUT_DIR."


```

Run the script
[Access eigenvec_to_csv.sh](/Wrappers/eigenvec_to_csv.sh)

```bash
sbash eigenvec_to_csv.sh
```

#### PCA

Move to the SCRIPTS directory

```bash
cd /scratch/MOryzae/SCRIPTS
```

Open nano text editor

```bash
nano pca_plot.sh
```

save the following sbatch script

```bash
#!/bin/bash

############# SLURM Configuration ##############

### Define Job name
#SBATCH --job-name=genome_pca_plot

### Define partition to use
#SBATCH -p normal

### Define number of CPUs to use
#SBATCH -c 8

### Specify the node to run on
#SBATCH --nodelist=node20

#################################################

########### Execution Command ###################

module load python/3.12.0  # Charge Python 3.12 sur le cluster

# Define directories
PCA_RESULTS_DIR="/scratch/MOryzae/PLINK"
OUTPUT_PLOT_DIR="/scratch/MOryzae/PLINK"

# List of directories to process
DIRECTORIES=("plink")

# Loop through each directory
for DIRECTORY in "${DIRECTORIES[@]}"; do
    PCA_FILE="$PCA_RESULTS_DIR/$DIRECTORY/dataset.eigenvec"
    OUTPUT_DIR="$OUTPUT_PLOT_DIR/$DIRECTORY"
    OUTPUT_PLOT_2D="$OUTPUT_DIR/dataset_2D.png"
    OUTPUT_PLOT_3D="$OUTPUT_DIR/dataset_3D.png"

    # Ensure the output directory exists
    mkdir -p "$OUTPUT_DIR"

    # Check if the PCA results file exists
    if [ -f "$PCA_FILE" ]; then
        echo "Processing PCA results for $DIRECTORY..."
        
        # Call the Python script to generate the plots
        python3 <<EOF
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

# Define input and output paths
pca_results_file = "$PCA_FILE"
output_plot_2D = "$OUTPUT_PLOT_2D"
output_plot_3D = "$OUTPUT_PLOT_3D"

# Read PCA results
try:
    pca_results = pd.read_csv(pca_results_file, sep=r'\s+', header=None)
    pca_results.columns = ['FID', 'IID', 'PC1', 'PC2', 'PC3']
except Exception as e:
    print(f"Error reading PCA results file {pca_results_file}: {e}")
    exit(1)

# Plot 2D scatter plot for PC1 vs PC2
try:
    plt.figure(figsize=(8, 6))
    plt.scatter(pca_results['PC1'], pca_results['PC2'], s=100)
    plt.title('PCA Results: $DIRECTORY (2D)')
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.grid()
    plt.savefig(output_plot_2D)
    plt.close()
except Exception as e:
    print(f"Error creating 2D plot for {pca_results_file}: {e}")
    exit(1)

# Plot 3D scatter plot for PC1, PC2, and PC3
try:
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    scatter = ax.scatter(pca_results['PC1'], pca_results['PC2'], pca_results['PC3'], s=100, c='blue', alpha=0.7)
    ax.set_title('PCA Results: $DIRECTORY (3D)')
    ax.set_xlabel('Principal Component 1')
    ax.set_ylabel('Principal Component 2')
    ax.set_zlabel('Principal Component 3')
    plt.savefig(output_plot_3D)
    plt.close()
except Exception as e:
    print(f"Error creating 3D plot for {pca_results_file}: {e}")
    exit(1)

# Verify plots were created
if not os.path.exists(output_plot_2D) or not os.path.exists(output_plot_3D):
    print(f"Error: Output plots not created for {pca_results_file}")
    exit(1)
EOF

        # Check if the Python script executed successfully
        if [ $? -eq 0 ]; then
            echo "Plots successfully created for $DIRECTORY: $OUTPUT_PLOT_2D, $OUTPUT_PLOT_3D"
        else
            echo "Error: Failed to create plots for $DIRECTORY."
            exit 1
        fi
    else
        echo "Error: PCA results file not found in $DIRECTORY."
        exit 1
    fi
done

echo "All PCA plots created successfully."

```

Run the script
[Access pca_plot.sh](/Wrappers/pca_plot.sh)

```bash
sbash pca_plot.sh
```

Using tools to interpret PCA results or plots

To do this, we have three frequently used tools at our disposal
 * Use k-means to partition the data, assuming a number of clusters.
 * Apply DBSCAN to detect dense clusters and identify outliers.
 * Use the elbow method and silhouette index to determine the optimal number of clusters.

Move to the SCRIPTS directory

```bash
cd /scratch/MOryzae/SCRIPTS
```

Open nano text editor

```bash
nano clustering_analysis.py
```

save the following sbatch script

```py
import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.cluster import KMeans, DBSCAN
from sklearn.metrics import silhouette_score

# Check for input arguments
if len(sys.argv) != 3:
    print("Usage: python clustering_analysis.py <PCA_FILE> <OUTPUT_DIR>")
    sys.exit(1)

# Input and output paths
pca_file = sys.argv[1]
output_dir = sys.argv[2]

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Load PCA data
try:
    print(f"Loading PCA data from {pca_file}...")
    data = pd.read_csv(pca_file, sep=r'\s+', header=None)
    data.columns = ['FID', 'IID', 'PC1', 'PC2', 'PC3']
    X = data[['PC1', 'PC2', 'PC3']]
except Exception as e:
    print(f"Error loading PCA data: {e}")
    sys.exit(1)

# Determine optimal number of clusters using the elbow method
print("Calculating optimal number of clusters (Elbow Method)...")
inertia = []
k_values = range(1, 10)
for k in k_values:
    kmeans = KMeans(n_clusters=k, random_state=42, n_init='auto')
    kmeans.fit(X)
    inertia.append(kmeans.inertia_)

plt.figure(figsize=(8, 5))
plt.plot(k_values, inertia, marker='o')
plt.title('Elbow Method')
plt.xlabel('Number of Clusters (k)')
plt.ylabel('Inertia')
elbow_plot_path = os.path.join(output_dir, 'elbow_method.png')
plt.savefig(elbow_plot_path)
print(f"Elbow Method plot saved to {elbow_plot_path}")

# Apply k-means clustering with k=3 (or adjust based on elbow results)
print("Applying k-means clustering...")
kmeans = KMeans(n_clusters=3, random_state=42, n_init='auto')
data['Cluster'] = kmeans.fit_predict(X)

# Sauvegarder les résultats avec les clusters dans un fichier CSV
output_csv = os.path.join(output_dir, 'clustered_data.csv')
data.to_csv(output_csv, index=False)
print(f"Données annotées avec clusters sauvegardées dans {output_csv}")


# Visualize k-means clusters in 3D
print("Generating 3D k-means plot...")
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
scatter = ax.scatter(data['PC1'], data['PC2'], data['PC3'], c=data['Cluster'], cmap='viridis', s=100)
plt.colorbar(scatter, ax=ax)
ax.set_title('k-means Clustering (k=3)')
ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
ax.set_zlabel('PC3')
kmeans_plot_path = os.path.join(output_dir, 'kmeans_pca_plot_3d.png')
plt.savefig(kmeans_plot_path)
print(f"3D k-means plot saved to {kmeans_plot_path}")

# Apply DBSCAN clustering
print("Applying DBSCAN clustering...")
dbscan = DBSCAN(eps=0.1, min_samples=3)
data['DBSCAN_Cluster'] = dbscan.fit_predict(X)

# Visualize DBSCAN clusters in 3D
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
scatter = ax.scatter(data['PC1'], data['PC2'], data['PC3'], c=data['DBSCAN_Cluster'], cmap='plasma', s=100)
plt.colorbar(scatter, ax=ax)
ax.set_title('DBSCAN Clustering')
ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
ax.set_zlabel('PC3')
dbscan_plot_path = os.path.join(output_dir, 'dbscan_pca_plot_3d.png')
plt.savefig(dbscan_plot_path)
print(f"3D DBSCAN plot saved to {dbscan_plot_path}")

# Calculate silhouette scores for different k
print("Calculating silhouette scores...")
silhouette_scores = []
for k in range(2, 10):
    kmeans = KMeans(n_clusters=k, random_state=42, n_init='auto')
    labels = kmeans.fit_predict(X)
    silhouette_scores.append(silhouette_score(X, labels))

plt.figure(figsize=(8, 5))
plt.plot(range(2, 10), silhouette_scores, marker='o')
plt.title('Silhouette Scores')
plt.xlabel('Number of Clusters (k)')
plt.ylabel('Silhouette Score')
silhouette_plot_path = os.path.join(output_dir, 'silhouette_scores.png')
plt.savefig(silhouette_plot_path)
print(f"Silhouette Scores plot saved to {silhouette_plot_path}")

print("Clustering analysis completed.")

```

Then

Open nano text editor

```bash
nano pca_learning.sh
```

save the following sbatch script

```bash
#!/bin/bash

############# SLURM Configuration ##############
#SBATCH --job-name=pca_learning
#SBATCH -p normal
#SBATCH -c 8
#SBATCH --nodelist=node20

#################################################

# Load necessary modules
module load python/3.12.0

# Define directories
PCA_RESULTS_DIR="/scratch/MOryzae/PLINK"
OUTPUT_PLOT_DIR="/scratch/MOryzae/PLINK/PCA"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_PLOT_DIR"

# List of subdirectories to process
DIRECTORIES=("plink")

# Loop over each directory
for DIRECTORY in "${DIRECTORIES[@]}"; do
    PCA_FILE="$PCA_RESULTS_DIR/$DIRECTORY/dataset.eigenvec"
    OUTPUT_DIR="$OUTPUT_PLOT_DIR"

    # Check if the PCA file exists
    if [[ -f "$PCA_FILE" ]]; then
        mkdir -p "$OUTPUT_DIR"
        echo "Processing PCA file: $PCA_FILE"

        # Run the Python script for clustering and plotting
        python3 clustering_analysis.py "$PCA_FILE" "$OUTPUT_DIR"

        echo "Processing completed for directory: $DIRECTORY"
    else
        echo "PCA file not found: $PCA_FILE"
    fi
done

```

Run the script
[Access pca_learning.sh](/Wrappers/pca_learning.sh)

```bash
sbash pca_learning.sh
```

Use generated files to interpret this step and plannig next analysis

For the next step,
Move the entire contents of the PLINK and PCA directories to the NAS

```bash
scp -r /scratch/MOryzae/PLINK/ san:/projects/medium/CIBiG_MOryzae/

```

Retrieve PLINK and PCA directories from the NAS on your local machine to analyze the results

```bash
scp -r login@bioinfo-san.ird.fr:/projects/medium/CIBiG_MOryzae/PLINK /path/to/working dorectory/on your laptop/

```

#### DAPC

Run the following commands on you local laptop

Open nano text editor

```bash
nano dapc_analysis.R
```

save the following sbatch script

```R
# Charger les bibliothèques nécessaires
if (!requireNamespace("adegenet")) install.packages("adegenet")
if (!requireNamespace("factoextra")) install.packages("factoextra") # Pour le clustering
if (!requireNamespace("ggplot2")) install.packages("ggplot2")
if (!requireNamespace("dplyr")) install.packages("dplyr")

library(adegenet)
library(factoextra)
library(ggplot2)
library(dplyr)

# Définir les chemins d'entrée et de sortie
input_file <- "/home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/PLINK/plink/dataset.eigenvec"
output_dir <- "/home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/PLINK/DAPC"

# Créer le répertoire de sortie si nécessaire
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Répertoire de sortie créé :", output_dir, "\n")
}

# Charger les données PCA
cat("Chargement des données PCA...\n")
pca_data <- read.table(input_file, header = FALSE)
colnames(pca_data) <- c("FID", "IID", "PC1", "PC2", "PC3")  # Modifier selon vos colonnes
X <- pca_data %>% select(PC1, PC2, PC3)

# Étape 1 : Clustering des individus (k-means)
cat("Application de k-means pour générer des groupes...\n")
set.seed(42)  # Pour assurer la reproductibilité
n_clusters <- 3  # Ajustez selon vos besoins ou utilisez la méthode du coude (voir ci-dessous)
kmeans_result <- kmeans(X, centers = n_clusters)

# Ajouter les groupes au jeu de données
pca_data$Group <- as.factor(kmeans_result$cluster)

# Sauvegarder les groupes dans un fichier CSV
groups_file <- file.path(output_dir, "groups_kmeans.csv")
write.csv(pca_data, groups_file, row.names = FALSE)
cat("Groupes sauvegardés dans :", groups_file, "\n")

# Étape optionnelle : Déterminer le nombre optimal de clusters
# Méthode du coude
cat("Déterminer le nombre optimal de clusters avec la méthode du coude...\n")
fviz_nbclust(X, kmeans, method = "wss") +
  labs(title = "Méthode du coude pour déterminer k")

# Étape 2 : DAPC
cat("Optimisation du nombre de PCs pour DAPC...\n")
dapc_initial <- dapc(X, pca_data$Group)
optimal_pcs <- optim.a.score(dapc_initial)
n_pcs <- optimal_pcs$n.pca
cat(paste("Nombre optimal de PCs :", n_pcs, "\n"))

# Réaliser la DAPC avec le nombre optimal de PCs
dapc_result <- dapc(X, pca_data$Group, n.pca = n_pcs)

# Étape 3 : Visualisation des résultats
cat("Génération du graphique des clusters DAPC...\n")
scatter_file <- file.path(output_dir, "dapc_scatter.png")
png(scatter_file, width = 800, height = 600)
scatter(dapc_result, scree.da = TRUE, posi.da = "bottomleft", scree.pca = TRUE)
dev.off()
cat("Graphique DAPC sauvegardé dans :", scatter_file, "\n")

# Graphique ggplot des clusters
dapc_df <- data.frame(dapc_result$ind.coord) %>%
  mutate(Group = pca_data$Group)  # Ajouter les groupes au dataframe

ggplot_file <- file.path(output_dir, "dapc_ggplot.png")
gg <- ggplot(dapc_df, aes(x = LD1, y = LD2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(title = "Clusters DAPC", x = "Discriminant Axis 1", y = "Discriminant Axis 2")
ggsave(ggplot_file, plot = gg, width = 8, height = 6)
cat("Graphique ggplot DAPC sauvegardé dans :", ggplot_file, "\n")

# Étape 4 : Contributions des variables
cat("Visualisation des contributions des variables...\n")
loading_file <- file.path(output_dir, "dapc_loadings.png")
png(loading_file, width = 800, height = 600)
loadingplot(dapc_result$var.contr, axis = 1, threshold = 0.005, lab.jitter = 1)
dev.off()
cat("Graphique des contributions sauvegardé dans :", loading_file, "\n")

# Étape 5 : Sauvegarder les résultats
results_file <- file.path(output_dir, "dapc_results_with_clusters.csv")
write.csv(dapc_df, results_file, row.names = FALSE)
cat("Résultats DAPC sauvegardés dans :", results_file, "\n")


```

Run the script
[Access dapc_analysis.R](/Wrappers/dapc_analysis.R)

```R
source("/path/to/working dorectory/on your laptop/dapc_analysis.R")
```

## 5. SNP FILTRATION

### 5.1. Evaluate some statistiques on VCF file

Move to the SCRIPTS directory

```bash
cd /scratch/MOryzae/SCRIPTS
```

Open nano text editor

```bash
nano snp_statistic.sh
```

save the following sbatch script

```bash
#!/bin/bash

############ SLURM Configuration ##############

### Define Job name
#SBATCH --job-name=snp_statistics

### Define partition to use
#SBATCH -p normal

### Define number of CPUs to use
#SBATCH -c 8

### Specify the node to run on
#SBATCH --nodelist=node20  # Run the job on node20

#################################################

########### Execution Commands ###################

# Define directories
SUBSET_VCF="/scratch/MOryzae/SNP/vcf_files/snp_correct.vcf.gz"
OUT_DIR="/scratch/MOryzae/SNP/Others_stats"

# Create output directory if it doesn't exist
mkdir -p "$OUT_DIR"

# Load necessary modules
module load vcftools/0.1.16

# Define base output name
OUT_PREFIX="${OUT_DIR}/output"

# Calculate allele frequency
vcftools --gzvcf "$SUBSET_VCF" --remove-indels --freq2 --max-alleles 2 --out "$OUT_PREFIX"

# Calculate mean depth per individual
vcftools --gzvcf "$SUBSET_VCF" --remove-indels --depth --out "$OUT_PREFIX"

# Calculate mean depth per site
vcftools --gzvcf "$SUBSET_VCF" --remove-indels --site-mean-depth --out "$OUT_PREFIX"

# Calculate site quality
vcftools --gzvcf "$SUBSET_VCF" --remove-indels --site-quality --out "$OUT_PREFIX"

# Calculate proportion of missing data per individual
vcftools --gzvcf "$SUBSET_VCF" --remove-indels --missing-indv --out "$OUT_PREFIX"

# Calculate proportion of missing data per site
vcftools --gzvcf "$SUBSET_VCF" --remove-indels --missing-site --out "$OUT_PREFIX"

# Calculate heterozygosity and inbreeding coefficient per individual
vcftools --gzvcf "$SUBSET_VCF" --remove-indels --het --out "$OUT_PREFIX"

echo "All SNP statistics have been computed. Results saved to $OUT_DIR"

```

Run the script
[Access snp_statistic.sh](/Wrappers/snp_statistic.sh)

```bash
sbash snp_statistic.sh
```

Use generated files to interpret this step and plannig next analysis

### 5.2. Examining statistics in R

Open nano text editor

```bash
nano snp_statistic.R
```

save the following sbatch script

```R
# Load necessary libraries
library(tidyverse)
library(ggplot2)

# Define working directories
input_dir <- "/home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/SNP/Others_stats"
output_dir <- "/home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/SNP/Others_stats/Plots"
summary_file <- file.path(output_dir, "summary_statistics.txt")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Set working directory
setwd(input_dir)

# Open summary file for writing
summary_conn <- file(summary_file, open = "w")

# Variant-based statistics

# Variant quality
var_qual <- read_delim(file.path(input_dir, "output.lqual"), delim = "\t",
                       col_names = c("chr", "pos", "qual"), skip = 1)

p <- ggplot(var_qual, aes(qual)) +
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light() +
  ggtitle("Variant Quality Distribution")
ggsave(file.path(output_dir, "variant_quality_density.png"), plot = p)

summary_var_qual <- summary(var_qual$qual)
writeLines("### Variant Quality Summary ###\n", summary_conn)
writeLines(capture.output(summary_var_qual), summary_conn)

# Variant mean depth
var_depth <- read_delim(file.path(input_dir, "output.ldepth.mean"), delim = "\t",
                        col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)

p <- ggplot(var_depth, aes(mean_depth)) +
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light() +
  xlim(0, 100) +
  ggtitle("Mean Depth per Variant")
ggsave(file.path(output_dir, "variant_depth_density.png"), plot = p)

summary_var_depth <- summary(var_depth$mean_depth)
writeLines("\n### Variant Mean Depth Summary ###\n", summary_conn)
writeLines(capture.output(summary_var_depth), summary_conn)

# Variant missingness
var_miss <- read_delim(file.path(input_dir, "output.lmiss"), delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)

p <- ggplot(var_miss, aes(fmiss)) +
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light() +
  ggtitle("Missingness per Variant")
ggsave(file.path(output_dir, "variant_missingness_density.png"), plot = p)

summary_var_miss <- summary(var_miss$fmiss)
writeLines("\n### Variant Missingness Summary ###\n", summary_conn)
writeLines(capture.output(summary_var_miss), summary_conn)

# Minor allele frequency
var_freq <- read_delim(file.path(input_dir, "output.frq"), delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)

# Calculate minor allele frequency
var_freq <- var_freq %>%
  mutate(maf = pmin(as.numeric(a1), as.numeric(a2), na.rm = TRUE))

p <- ggplot(var_freq, aes(maf)) +
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light() +
  ggtitle("Minor Allele Frequency Distribution")
ggsave(file.path(output_dir, "maf_density.png"), plot = p)

summary_var_freq <- summary(var_freq$maf)
writeLines("\n### Minor Allele Frequency Summary ###\n", summary_conn)
writeLines(capture.output(summary_var_freq), summary_conn)

# Individual-based statistics

# Mean depth per individual
ind_depth <- read_delim(file.path(input_dir, "output.idepth"), delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)

p <- ggplot(ind_depth, aes(depth)) +
  geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light() +
  ggtitle("Mean Depth per Individual")
ggsave(file.path(output_dir, "individual_depth_histogram.png"), plot = p)

summary_ind_depth <- summary(ind_depth$depth)
writeLines("\n### Individual Mean Depth Summary ###\n", summary_conn)
writeLines(capture.output(summary_ind_depth), summary_conn)

# Missing data per individual
ind_miss <- read_delim(file.path(input_dir, "output.imiss"), delim = "\t",
                       col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)

p <- ggplot(ind_miss, aes(fmiss)) +
  geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light() +
  ggtitle("Missing Data per Individual")
ggsave(file.path(output_dir, "individual_missing_data_histogram.png"), plot = p)

summary_ind_miss <- summary(ind_miss$fmiss)
writeLines("\n### Individual Missing Data Summary ###\n", summary_conn)
writeLines(capture.output(summary_ind_miss), summary_conn)

# Heterozygosity and inbreeding coefficient
ind_het <- read_delim(file.path(input_dir, "output.het"), delim = "\t",
                      col_names = c("ind", "ho", "he", "nsites", "f"), skip = 1)

p <- ggplot(ind_het, aes(f)) +
  geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light() +
  ggtitle("Inbreeding Coefficient Distribution")
ggsave(file.path(output_dir, "inbreeding_coefficient_histogram.png"), plot = p)

summary_ind_het <- summary(ind_het$f)
writeLines("\n### Inbreeding Coefficient Summary ###\n", summary_conn)
writeLines(capture.output(summary_ind_het), summary_conn)

# Close summary file
close(summary_conn)

message("All plots and summary statistics have been saved to ", output_dir)

```

Run the script
[Access snp_statistic.R](/Wrappers/snp_statistic.R)

```bash
source("/path/to/working dorectory/on your laptop/snp_statistic.R")
```

Use generated files to interpret this step and plannig next analysis


### 5.3. String mode SNP filtering

Move to the SCRIPTS directory

```bash
cd /scratch/MOryzae/SCRIPTS
```

Open nano text editor

```bash
nano snp_filtering.sh
```

save the following sbatch script

```bash
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

```

Run the script
[Access snp_filtering.sh](/Wrappers/snp_filtering.sh)

```bash
sbash snp_filtering.sh
```

Use generated files to interpret this step and plannig next analysis

### 5.4. PCA & DAPC

#### Generate PCA using genotyping information contained in VCF

Plink alllows to create a PCA (principal components analysis) of samples, so that we can easily evaluate genetic distance between samples.

This will generate a matrix of coordinates in the different component. By default, it provides the first 20 principal components of the variance-standardized relationship matrix. We will focus only the first 3 axes for subsequent visualization (--pca 3)

Create a PLINK directory in the working directory

```bash
mkdir -p /scratch/MOryzae/PLINK2
```

Move to the created directory

```bash
cd PLINK2/
```

Load plink module first

```bash
module load plink/1.9
```

Start evaluation of missing data

```bash
plink -vcf /scratch/MOryzae/SNP/vcf_filtered/filtered_snps.vcf.gz --allow-extra-chr --cluster --matrix --pca 3 --mind --out ./plink/dataset
```


#### Convert "eigenvec" generate file to csv file

Move to the SCRIPTS directory

```bash
cd /scratch/MOryzae/SCRIPTS
```

Open nano text editor

```bash
nano eigenvec2_to_csv.sh
```

save the following sbatch script

```bash
#!/bin/bash

############# SLURM Configuration ##############

### Define Job name
#SBATCH --job-name=eigenvec_to_csv

### Define partition to use
#SBATCH -p normal

### Define number of CPUs to use
#SBATCH -c 8

### Specify the node to run on
#SBATCH --nodelist=node20  # Run the job on node20

#################################################

########### Execution Command ###################

# Define directories
INPUT_DIR="/scratch/MOryzae/PLINK2"
OUTPUT_DIR="/scratch/MOryzae/PLINK2"

# List of directories to process
DIRECTORIES=("plink")

# Loop through each directory
for DIR in "${DIRECTORIES[@]}"; do
    PCA_FILE="$INPUT_DIR/$DIR/dataset.eigenvec"
    OUTPUT_CSV="$OUTPUT_DIR/$DIR/dataset.csv"

    # Check if PCA results exist
    if [ -f "$PCA_FILE" ]; then
        echo "Processing PCA results for $DIR..."

        # Convert eigenvec file to CSV format
        awk 'NR==1{print "FID,IID,PC1,PC2,PC3"} NR>1{print $1","$2","$3","$4","$5}' "$PCA_FILE" > "$OUTPUT_CSV"

        # Check if conversion was successful
        if [ $? -eq 0 ]; then
            echo "Conversion successful: $OUTPUT_CSV created."
        else
            echo "Error: Failed to convert $PCA_FILE to CSV."
            exit 1
        fi
    else
        echo "Error: PCA results file not found in $DIR."
        exit 1
    fi
done

echo "PCA analysis completed. Results are saved in $OUTPUT_DIR."


```

Run the script
[Access eigenvec2_to_csv.sh](/Wrappers/eigenvec2_to_csv.sh)

```bash
sbash eigenvec2_to_csv.sh
```

#### PCA

Move to the SCRIPTS directory

```bash
cd /scratch/MOryzae/SCRIPTS
```

Open nano text editor

```bash
nano pca2_plot.sh
```

save the following sbatch script

```bash
#!/bin/bash

############# SLURM Configuration ##############

### Define Job name
#SBATCH --job-name=genome_pca_plot

### Define partition to use
#SBATCH -p normal

### Define number of CPUs to use
#SBATCH -c 8

### Specify the node to run on
#SBATCH --nodelist=node20

#################################################

########### Execution Command ###################

module load python/3.12.0  # Charge Python 3.12 sur le cluster

# Define directories
PCA_RESULTS_DIR="/scratch/MOryzae/PLINK2"
OUTPUT_PLOT_DIR="/scratch/MOryzae/PLINK2"

# List of directories to process
DIRECTORIES=("plink")

# Loop through each directory
for DIRECTORY in "${DIRECTORIES[@]}"; do
    PCA_FILE="$PCA_RESULTS_DIR/$DIRECTORY/dataset.eigenvec"
    OUTPUT_DIR="$OUTPUT_PLOT_DIR/$DIRECTORY"
    OUTPUT_PLOT_2D="$OUTPUT_DIR/dataset_2D.png"
    OUTPUT_PLOT_3D="$OUTPUT_DIR/dataset_3D.png"

    # Ensure the output directory exists
    mkdir -p "$OUTPUT_DIR"

    # Check if the PCA results file exists
    if [ -f "$PCA_FILE" ]; then
        echo "Processing PCA results for $DIRECTORY..."
        
        # Call the Python script to generate the plots
        python3 <<EOF
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

# Define input and output paths
pca_results_file = "$PCA_FILE"
output_plot_2D = "$OUTPUT_PLOT_2D"
output_plot_3D = "$OUTPUT_PLOT_3D"

# Read PCA results
try:
    pca_results = pd.read_csv(pca_results_file, sep=r'\s+', header=None)
    pca_results.columns = ['FID', 'IID', 'PC1', 'PC2', 'PC3']
except Exception as e:
    print(f"Error reading PCA results file {pca_results_file}: {e}")
    exit(1)

# Plot 2D scatter plot for PC1 vs PC2
try:
    plt.figure(figsize=(8, 6))
    plt.scatter(pca_results['PC1'], pca_results['PC2'], s=100)
    plt.title('PCA Results: $DIRECTORY (2D)')
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.grid()
    plt.savefig(output_plot_2D)
    plt.close()
except Exception as e:
    print(f"Error creating 2D plot for {pca_results_file}: {e}")
    exit(1)

# Plot 3D scatter plot for PC1, PC2, and PC3
try:
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    scatter = ax.scatter(pca_results['PC1'], pca_results['PC2'], pca_results['PC3'], s=100, c='blue', alpha=0.7)
    ax.set_title('PCA Results: $DIRECTORY (3D)')
    ax.set_xlabel('Principal Component 1')
    ax.set_ylabel('Principal Component 2')
    ax.set_zlabel('Principal Component 3')
    plt.savefig(output_plot_3D)
    plt.close()
except Exception as e:
    print(f"Error creating 3D plot for {pca_results_file}: {e}")
    exit(1)

# Verify plots were created
if not os.path.exists(output_plot_2D) or not os.path.exists(output_plot_3D):
    print(f"Error: Output plots not created for {pca_results_file}")
    exit(1)
EOF

        # Check if the Python script executed successfully
        if [ $? -eq 0 ]; then
            echo "Plots successfully created for $DIRECTORY: $OUTPUT_PLOT_2D, $OUTPUT_PLOT_3D"
        else
            echo "Error: Failed to create plots for $DIRECTORY."
            exit 1
        fi
    else
        echo "Error: PCA results file not found in $DIRECTORY."
        exit 1
    fi
done

echo "All PCA plots created successfully."

```

Run the script
[Access pca2_plot.sh](/Wrappers/pca2_plot.sh)

```bash
sbash pca2_plot.sh
```

Using tools to interpret PCA results or plots

To do this, we have three frequently used tools at our disposal
 * Use k-means to partition the data, assuming a number of clusters.
 * Apply DBSCAN to detect dense clusters and identify outliers.
 * Use the elbow method and silhouette index to determine the optimal number of clusters.

Open nano text editor

```bash
nano pca2_learning.sh
```

save the following sbatch script

```bash
#!/bin/bash

############# SLURM Configuration ##############
#SBATCH --job-name=pca_learning
#SBATCH -p normal
#SBATCH -c 8
#SBATCH --nodelist=node20

#################################################

# Load necessary modules
module load python/3.12.0

# Define directories
PCA_RESULTS_DIR="/scratch/MOryzae/PLINK2"
OUTPUT_PLOT_DIR="/scratch/MOryzae/PLINK2/PCA"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_PLOT_DIR"

# List of subdirectories to process
DIRECTORIES=("plink")

# Loop over each directory
for DIRECTORY in "${DIRECTORIES[@]}"; do
    PCA_FILE="$PCA_RESULTS_DIR/$DIRECTORY/dataset.eigenvec"
    OUTPUT_DIR="$OUTPUT_PLOT_DIR"

    # Check if the PCA file exists
    if [[ -f "$PCA_FILE" ]]; then
        mkdir -p "$OUTPUT_DIR"
        echo "Processing PCA file: $PCA_FILE"

        # Run the Python script for clustering and plotting
        python3 clustering_analysis.py "$PCA_FILE" "$OUTPUT_DIR"

        echo "Processing completed for directory: $DIRECTORY"
    else
        echo "PCA file not found: $PCA_FILE"
    fi
done

```

Run the script
[Access pca2_learning.sh](/Wrappers/pca2_learning.sh)

```bash
sbash pca2_learning.sh
```

Use generated files to interpret this step and plannig next analysis

For the next step,
Move the entire contents of the PLINK and PCA directories to the NAS

```bash
scp -r /scratch/MOryzae/PLINK2/ san:/projects/medium/CIBiG_MOryzae/

```

Retrieve PLINK and PCA directories from the NAS on your local machine to analyze the results

```bash
scp -r login@bioinfo-san.ird.fr:/projects/medium/CIBiG_MOryzae/PLINK2 /path/to/working dorectory/on your laptop/

```

#### DAPC

Run the following commands on you local laptop

Open nano text editor

```bash
nano dapc2_analysis.R
```

save the following sbatch script

```R
# Charger les bibliothèques nécessaires
if (!requireNamespace("adegenet")) install.packages("adegenet")
if (!requireNamespace("factoextra")) install.packages("factoextra") # Pour le clustering
if (!requireNamespace("ggplot2")) install.packages("ggplot2")
if (!requireNamespace("dplyr")) install.packages("dplyr")

library(adegenet)
library(factoextra)
library(ggplot2)
library(dplyr)

# Définir les chemins d'entrée et de sortie
input_file <- "/home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/PLINK2/plink/dataset.eigenvec"
output_dir <- "/home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/PLINK2/DAPC"

# Créer le répertoire de sortie si nécessaire
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Répertoire de sortie créé :", output_dir, "\n")
}

# Charger les données PCA
cat("Chargement des données PCA...\n")
pca_data <- read.table(input_file, header = FALSE)
colnames(pca_data) <- c("FID", "IID", "PC1", "PC2", "PC3")  # Modifier selon vos colonnes
X <- pca_data %>% select(PC1, PC2, PC3)

# Étape 1 : Clustering des individus (k-means)
cat("Application de k-means pour générer des groupes...\n")
set.seed(42)  # Pour assurer la reproductibilité
n_clusters <- 3  # Ajustez selon vos besoins ou utilisez la méthode du coude (voir ci-dessous)
kmeans_result <- kmeans(X, centers = n_clusters)

# Ajouter les groupes au jeu de données
pca_data$Group <- as.factor(kmeans_result$cluster)

# Sauvegarder les groupes dans un fichier CSV
groups_file <- file.path(output_dir, "groups_kmeans.csv")
write.csv(pca_data, groups_file, row.names = FALSE)
cat("Groupes sauvegardés dans :", groups_file, "\n")

# Étape optionnelle : Déterminer le nombre optimal de clusters
# Méthode du coude
cat("Déterminer le nombre optimal de clusters avec la méthode du coude...\n")
fviz_nbclust(X, kmeans, method = "wss") +
  labs(title = "Méthode du coude pour déterminer k")

# Étape 2 : DAPC
cat("Optimisation du nombre de PCs pour DAPC...\n")
dapc_initial <- dapc(X, pca_data$Group)
optimal_pcs <- optim.a.score(dapc_initial)
n_pcs <- optimal_pcs$n.pca
cat(paste("Nombre optimal de PCs :", n_pcs, "\n"))

# Réaliser la DAPC avec le nombre optimal de PCs
dapc_result <- dapc(X, pca_data$Group, n.pca = n_pcs)

# Étape 3 : Visualisation des résultats
cat("Génération du graphique des clusters DAPC...\n")
scatter_file <- file.path(output_dir, "dapc_scatter.png")
png(scatter_file, width = 800, height = 600)
scatter(dapc_result, scree.da = TRUE, posi.da = "bottomleft", scree.pca = TRUE)
dev.off()
cat("Graphique DAPC sauvegardé dans :", scatter_file, "\n")

# Graphique ggplot des clusters
dapc_df <- data.frame(dapc_result$ind.coord) %>%
  mutate(Group = pca_data$Group)  # Ajouter les groupes au dataframe

ggplot_file <- file.path(output_dir, "dapc_ggplot.png")
gg <- ggplot(dapc_df, aes(x = LD1, y = LD2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(title = "Clusters DAPC", x = "Discriminant Axis 1", y = "Discriminant Axis 2")
ggsave(ggplot_file, plot = gg, width = 8, height = 6)
cat("Graphique ggplot DAPC sauvegardé dans :", ggplot_file, "\n")

# Étape 4 : Contributions des variables
cat("Visualisation des contributions des variables...\n")
loading_file <- file.path(output_dir, "dapc_loadings.png")
png(loading_file, width = 800, height = 600)
loadingplot(dapc_result$var.contr, axis = 1, threshold = 0.005, lab.jitter = 1)
dev.off()
cat("Graphique des contributions sauvegardé dans :", loading_file, "\n")

# Étape 5 : Sauvegarder les résultats
results_file <- file.path(output_dir, "dapc_results_with_clusters.csv")
write.csv(dapc_df, results_file, row.names = FALSE)
cat("Résultats DAPC sauvegardés dans :", results_file, "\n")


```

Run the script
[Access dapc2_analysis.R](/Wrappers/dapc2_analysis.R)

```R
source("/path/to/working dorectory/on your laptop/dapc2_analysis.R")
```
