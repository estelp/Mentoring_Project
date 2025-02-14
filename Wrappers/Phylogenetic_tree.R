library(vcfR)
library(adegenet)
library(poppr)
library(ape)
library(RColorBrewer)
library(dartR)

# Chargement du fichier VCF
vcf_file <- "/home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/SNP/vcf_filtered/filtered_snps.vcf"
rubi.VCF <- read.vcfR(vcf_file)

# Chargement des métadonnées de population avec séparation par virgule
pop_file <- "/home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/SNP/vcf_filtered/population_data.csv"
pop.data <- read.table(pop_file, header = TRUE, sep = ",", stringsAsFactors = FALSE)

# Vérification des premières lignes du fichier pour confirmer la structure
head(pop.data)

# Vérification des colonnes
if (!all(c("AccessID", "Population") %in% colnames(pop.data))) {
  stop("Le fichier population_data.csv ne contient pas les colonnes AccessID et Population.")
}

# Vérification des correspondances entre le VCF et les métadonnées
if (!all(colnames(rubi.VCF@gt)[-1] %in% pop.data$AccessID)) {
  stop("Certains échantillons du VCF ne sont pas présents dans le fichier de population.")
}

# Réorganisation des données pour correspondre au VCF
pop.data <- pop.data[match(colnames(rubi.VCF@gt)[-1], pop.data$AccessID), ]

# Conversion du VCF en objet genlight
gl.rubi <- vcfR2genlight(rubi.VCF)

# Définition de la ploïdie et des groupes de population
ploidy(gl.rubi) <- 1
pop(gl.rubi) <- pop.data$Population

# Calcul des distances génétiques
x.dist <- poppr::bitwise.dist(gl.rubi)

# Construction des arbres phylogénétiques avec bootstrap (1 000 réplicats)
set.seed(123)  # Assurer la reproductibilité
tree_upgma <- aboot(gl.rubi, tree = "upgma", distance = bitwise.dist, sample = 500, 
                    showtree = FALSE, cutoff = 50, quiet = TRUE)

tree_nj <- aboot(gl.rubi, tree = "nj", distance = bitwise.dist, sample = 500, 
                 showtree = FALSE, cutoff = 50, quiet = TRUE)

# Définition des couleurs pour les populations
cols <- brewer.pal(n = length(unique(pop(gl.rubi))), name = "Dark2")

# Affichage de l'arbre UPGMA
plot.phylo(tree_upgma, cex = 0.8, font = 2, adj = 0, tip.color = cols[as.factor(pop(gl.rubi))], 
           main = "Arbre UPGMA basé sur les SNPs")

# Affichage de l'arbre NJ
plot.phylo(tree_nj, cex = 0.8, font = 2, adj = 0, tip.color = cols[as.factor(pop(gl.rubi))], 
           main = "Arbre NJ basé sur les SNPs")

# Sauvegarde des arbres au format Newick
write.tree(tree_upgma, file = "/home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/SNP/tree/tree_upgma.newick")
write.tree(tree_nj, file = "/home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/SNP/tree/tree_nj.newick")


