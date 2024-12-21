# Installer et charger les packages nécessaires
install.packages("ggplot2")
install.packages("sf")  # Pour manipuler les données géographiques
install.packages("rnaturalearth")  # Pour obtenir les données géographiques mondiales

library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

# Télécharger les données géographiques des pays
world <- ne_countries(scale = "medium", returnclass = "sf")

# Lister les pays
pays_liste <- world$name

# Afficher la liste des pays
print(pays_liste)

# Filtrer les pays d'intérêt
countries_of_interest <- c("Côte d'Ivoire", "Tanzania", "Cameroon", "Mali", "Iran", "Japan", 
                           "United States of America", "India", "Bangladesh", "China", "Bolivia", 
                           "Brazil", "Paraguay", "Nigeria", "Togo", "Gabon", "Nepal", "Vietnam", 
                           "Benin", "Argentina", "France", "Burkina Faso", "Philippines", "Thailand", 
                           "Turkey", "Italy")

# Filtrer le monde pour les pays sélectionnés
world_filtered <- world[world$name %in% countries_of_interest, ]

# Générer une couleur aléatoire pour chaque pays d'intérêt
set.seed(123)  # Fixer la graine pour la reproductibilité
world_filtered$color <- sample(colors(), nrow(world_filtered), replace = TRUE)

# Ajouter une couleur neutre (par exemple gris clair) aux autres pays
world$color <- ifelse(world$name %in% countries_of_interest, 
                      world_filtered$color,  # couleur aléatoire pour les pays d'intérêt
                      "lightgrey")           # couleur neutre pour les autres pays

# Créer la carte avec un arrière-plan vide
ggplot(data = world) +
  geom_sf(aes(fill = color), color = "black") +  # Ajouter la couche de base avec tous les pays
  geom_sf(data = world_filtered, aes(fill = color), color = "black", size = 0.5) +  # Mettre en surbrillance les pays sélectionnés
  scale_fill_identity() +  # Utiliser les couleurs comme elles sont
  theme_void() +  # Appliquer un thème sans fond ni axes (arrière-plan vide)
  labs(title = "Maping of countries") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.position = "right")  # Positionner la légende à droite
