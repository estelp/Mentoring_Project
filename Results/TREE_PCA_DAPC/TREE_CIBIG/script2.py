import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

# Utiliser un backend non interactif pour éviter les erreurs d'affichage
matplotlib.use('Agg')  # Utilise le backend 'Agg' pour générer des fichiers image

# Charger les données depuis le fichier Excel
data = pd.read_excel("/media/name/ESTEL3/TREE_CIBIG/K_file.xlsx")

# Convertir les données des clusters en valeurs numériques (au cas où il y aurait des valeurs non numériques)
clusters = data.columns[1:]  # Ignorer la première colonne (les isolats)
data[clusters] = data[clusters].apply(pd.to_numeric, errors='coerce')  # Convertir en numériques

# Préparer les données pour le graphique
isolates = data["Isolate"]
num_clusters = len(clusters)

# Définir une palette de couleurs
colors = plt.cm.get_cmap("tab10", num_clusters)  # Définir une palette de couleurs

# Créer une figure et un axe pour le graphique
fig, ax = plt.subplots(figsize=(10, 6))

# Initialiser la base comme un tableau de zéros (numpy array)
bottoms = [0] * len(isolates)  # Cela était une liste, la conversion en numpy array résout l'erreur

# Créer les barres empilées
for i, cluster in enumerate(clusters):
    ax.bar(isolates, data[cluster], bottom=bottoms, label=cluster, color=colors(i))
    bottoms = [x + y for x, y in zip(bottoms, data[cluster])]  # Mettre à jour la base

# Personnaliser le graphique
plt.xticks(rotation=90)
plt.legend(title="Clusters")
plt.ylabel("Proportion")
plt.title("Clustering Results")
plt.tight_layout()

# Enregistrer le graphique dans un fichier PNG
plt.savefig("/media/name/ESTEL3/TREE_CIBIG/clustering_stacked_bar.png")

# Afficher un message pour indiquer que l'image a été enregistrée
print("Graphique enregistré dans : /media/name/ESTEL3/TREE_CIBIG/clustering_stacked_bar.png")

