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
* <u>FASTQC</u>

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
#SBATCH --nodelist=node13  # Spécifie que le job doit être exécuté sur node13

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
[Access FastqQC.sh]/Wrappers/FastQC)

```bash
sbash FastQC.sh
```

At the end of the task, check the contents

```bash
ls -lh /scratch/MOryzae/QC/FastQC/
```

### Use generated html files to check read quality


* <u>MULTIQC</u>

### 3. Mapping
### 4. SNP calling
