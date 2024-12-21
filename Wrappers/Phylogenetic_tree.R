library(vcfR)
library(ade4)
library(adegenet)
library(poppr)
library(ape)
library(RColorBrewer)
library(trekcolors)
library('dartR')

#Aller dans le repertoire
rubi.VCF <- read.vcfR("/home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/SNP/vcf_filtered/filtered_snps.vcf")
pop.data <- read.table("/home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/SNP/vcf_filtered/filtered_snps.vcf", sep = "\t", header = TRUE)
all(colnames(rubi.VCF@gt)[-1] == pop.data$AccessID)

#VCF en objet genlight.
gl.rubi <- vcfR2genlight(rubi.VCF)

#################################Apparition de waring
#genlight ne supporte que les loci bialléliques (avec deux allèles) alors que le fichier contient des loci multialléliques
#genlight omet les locis multialleliques 
################################# Pas d'effet

ploidy(gl.rubi) <- 1
pop(gl.rubi) <- pop.data$State
gl.rubi
x.dist <- dist(gl.rubi)
x.dist <- poppr::bitwise.dist(gl.rubi)
########################################################
# Charger les bibliothèques nécessaires
library(ape)
# Construire l'arbre avec la méthode UPGMA et la distance bitwise
tree <- aboot(gl.rubi, tree = "upgma", distance = bitwise.dist, sample = 100, 
              showtree = F, cutoff = 50, quiet = T)
cols <- brewer.pal(n = nPop(gl.rubi), name = "Dark2") ##Optionnel Itol regle ce warning 
# Tracer l'arbre avec des paramètres de mise en forme
plot.phylo(tree, cex = 0.8, font = 2, adj = 0, tip.color = cols[pop(gl.rubi)])
# Exporter l'arbre au format Newick dans un fichier
write.tree(tree, file = "tree_output.newick")
