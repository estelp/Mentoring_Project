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
input_file <- "/home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/PLINK/only_snp/dataset.eigenvec"
output_dir <- "/home/name/Documents/Projet_CIBiG/Mentoring_Project/Results/PLINK/DAPC/only_snp"

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

