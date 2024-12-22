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
##Option Itol avec newick regle la notation, la coloration ect.. 

# Création d'un nom de fichier plus descriptif
output_file <- file.path("/home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/TREE_PCA_DAPC/", 
                         paste0("tree_output_", Sys.Date(), ".newick"))

# Exportation de l'arbre en format Newick avec le chemin complet
write.tree(tree, file = output_file)

# Affichage d'un message de confirmation
cat("L'arbre a été exporté avec succès vers : ", output_file, "\n")

################################################################
#Principal components analysis
################################################################

#En fonction de l'abre phylogentique , Il faut definir Pop_data en creant la corespondance entre le nom de l'individu et le groupe de population auquel il appartient

pop_data <- c("Pop4","Pop6", "Pop6", "Pop6", "Pop6", "Pop6", "Pop1", "Pop5", "Pop1", "Pop3", 
              "Pop3", "Pop4", "Pop1", "Pop6", "Pop6", "Pop3", "Pop3", "Pop7", "Pop7","Pop3",
              "Pop7","Pop4","Pop7","Pop7", "Pop6","Pop6", "Pop1", "Pop6", "Pop6", "Pop6", 
              "Pop6", "Pop5", "Pop6", "Pop5", "Pop1", "Pop6", "Pop4", "Pop3", "Pop2", "Pop3", 
              "Pop3", "Pop3", "Pop3", "Pop3", "Pop3", "Pop3", "Pop4", "Pop4", "Pop5","Pop5",
              "Pop5", "Pop5", "Pop5", "Pop5", "Pop5", "Pop4", "Pop5", "Pop6", "Pop1","Pop3",
              "Pop1", "Pop1", "Pop1", "Pop7","Pop6","Pop6", "Pop6", "Pop6", "Pop6", "Pop6",
              "Pop6", "Pop6","Pop6", "Pop6", "Pop6", "Pop3", "Pop3","Pop2", "Pop2","Pop3",
              "Pop3", "Pop2", "Pop2", "Pop4", "Pop3", "Pop5", "Pop1", "Pop1","Pop6")



################################################################Good
# Vérifiez que le nombre d'éléments dans pop_data est bien égal au nombre d'individus
length(pop_data) == nInd(gl.rubi)  # Devrait renvoyer TRUE
# Calculer l'ACP à partir des données génétiques (gl.rubi) avec 3 composantes principales (PC1, PC2, PC3)
rubi.pca <- glPca(gl.rubi, nf = 7)

# Vérifier les valeurs propres (eig) et le pourcentage de variance expliqué
barplot(100 * rubi.pca$eig / sum(rubi.pca$eig), col = heat.colors(50), main = "PCA Eigenvalues")
title(ylab = "Percent of variance\nexplained", line = 2)
title(xlab = "Eigenvalues", line = 1)

# Conversion des scores PCA en data frame
rubi.pca.scores <- as.data.frame(rubi.pca$scores)

# Ajouter les informations de populations/clusters au score PCA
rubi.pca.scores$pop <- pop_data  # Si pop_data contient les informations des populations des individus

# Afficher les résultats de l'ACP : Première étape de PCA (PC1 et PC2) sur un graphique 2D
library(ggplot2)
set.seed(9)  # Assurer la reproductibilité des résultats

# Créer le graphique avec ggplot2 pour les deux premières composantes principales (PC1, PC2)
p <- ggplot(rubi.pca.scores, aes(x = PC1, y = PC2, color = pop)) +
  geom_point(size = 2) +  # Ajouter des points représentant chaque individu
  scale_color_manual(values = c("tomato", "lightblue", "olivedrab", "gold", "purple", "cyan", "orange")) +  # Définir les couleurs des populations
  geom_hline(yintercept = 0) +  # Ajouter une ligne horizontale à y = 0
  geom_vline(xintercept = 0) +  # Ajouter une ligne verticale à x = 0
  theme_bw() +  # Appliquer un fond blanc avec des axes noirs
  labs(title = "PCA - First two components", x = "PC1", y = "PC2") # Ajouter un titre et des labels

# Affichage du graphique
print(p)

# Optionnel : Ajouter une ellipse autour des clusters (95% de confiance)
p <- p + stat_ellipse(level = 0.95, size = 0.5)
print(p)
##################################################################
# glPca effectue une PCA, nf = 3 : calculer les trois parametre (PC1, PC2, et PC3) 
#rubi.pca$eig : Cette ligne récupère les valeurs propres (eigenvalues) des composantes principales de votre PCA. 
#Les valeurs propres indiquent la quantité de variance expliquée par chaque composante principale.
#100*rubi.pca$eig/sum(rubi.pca$eig) : le pourcentage de variance expliqué par chaque composante principale
#Conversion des scores PCA en data frame et ajout des clusters
#la position de chaque observation (position d'individu) sur chaque composant principal
#pop, qui contient les informations de populations ou clusters des individus dans gl.rubi. 
#Ces clusters seront utilisés pour colorier les points dans le graphique et pour délimiter les groupes avec des ellipses
#set.seed(9)#permet de garantir la reproductibilité des résultats.
# La couleur des points est déterminée par la colonne pop, qui représente les clusters ou populations.
# Ajout des points sur le graphique
#stat_ellipse(level = 0.95) : trace des ellipses autour des groupes d'individus pour indiquer l'étendue de chaque cluster dans l'espace PCA. 
#L'argument level = 0.95 signifie que l'ellipse couvrira 95 % des données de chaque groupe.
#scale_color_manual(values = cols) : définir manuellement les couleurs à utiliser pour chaque groupe (populations ou clusters).
#La variable cols contient une palette de couleurs personnalisée, que vous devez définir avant d'exécuter ce code. Par exemple, 
#ajoute une ligne horizontale à y = 0.
#Ajoute une ligne verticale à x = 0, 
#theme_bw()#Applique un fond blanc et des axes noirs, pour améliorer la lisibilité.

#############################################################################
#Discriminant Principal components analysis
#############################################################################

library(adegenet)

# Si vous avez des données génétiques dans un objet 'genind', vous pouvez utiliser dapc() comme suit :
# Exemple de création d'un objet genind à partir de genlight
# Convertir l'objet genlight en une matrice
genlight_matrix <- as.matrix(gl.rubi)

# Assurez-vous que la matrice a bien des lignes et des colonnes définies
dimnames(genlight_matrix) <- list(indNames(gl.rubi), locNames(gl.rubi))

# Convertir la matrice en un objet genind
genind_rubi <- as.genind(genlight_matrix)

# Assigner les populations à l'objet genind
pop(genind_rubi) <- c("Pop4","Pop6", "Pop6", "Pop6", "Pop6", "Pop6", "Pop1", "Pop5", "Pop1", "Pop3", 
                      "Pop3", "Pop4", "Pop1", "Pop6", "Pop6", "Pop3", "Pop3", "Pop7", "Pop7","Pop3",
                      "Pop7","Pop4","Pop7","Pop7", "Pop6","Pop6", "Pop1", "Pop6", "Pop6", "Pop6", 
                      "Pop6", "Pop5", "Pop6", "Pop5", "Pop1", "Pop6", "Pop4", "Pop3", "Pop2", "Pop3", 
                      "Pop3", "Pop3", "Pop3", "Pop3", "Pop3", "Pop3", "Pop4", "Pop4", "Pop5","Pop5",
                      "Pop5", "Pop5", "Pop5", "Pop5", "Pop5", "Pop4", "Pop5", "Pop6", "Pop1","Pop3",
                      "Pop1", "Pop1", "Pop1", "Pop7","Pop6","Pop6", "Pop6", "Pop6", "Pop6", "Pop6",
                      "Pop6", "Pop6","Pop6", "Pop6", "Pop6", "Pop3", "Pop3","Pop2", "Pop2","Pop3",
                      "Pop3", "Pop2", "Pop2", "Pop4", "Pop3", "Pop5", "Pop1", "Pop1","Pop6")
# Effectuer l'analyse DAPC
pnw.dapc <- dapc(genind_rubi, n.pca = 7, n.da = 2)
            scale_color_manual(values = c("tomato", "lightblue", "olivedrab", "gold", "purple", "cyan", "orange"))

# Résumé des résultats
summary(pnw.dapc)

# Visualisation des résultats (scatter plot)
scatter(pnw.dapc)

# Visualiser la projection des individus sur les deux premiers axes
scatter(pnw.dapc, posi.da = "bottomleft", bg = "white")

#L'emplacement assign.per.pop indique les proportions de réaffectation réussie (sur la base des fonctions discriminantes) des individus vers leurs groupes d'origine
assignplot(pnw.dapc, subset=1:89)

###################################################################################
