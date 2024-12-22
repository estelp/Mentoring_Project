import pandas as pd
import matplotlib.pyplot as plt

# Chemin vers le fichier Excel
file_path = "/media/name/ESTEL3/TREE_CIBIG/K_file.xlsx"

# Charger les données depuis le fichier Excel
data = pd.read_excel(file_path)

# Préparer les données pour le graphique
isolates = data["Isolate"]
clusters = data.columns[1:]
colors = data.iloc[:, 1:]

# Création du graphique en barres empilées
fig, ax = plt.subplots(figsize=(12, 8))

# Boucle pour chaque cluster (colonnes après "Isolate") et ajout des barres empilées
for i, cluster in enumerate(clusters):
    bottom = colors.iloc[:, :i].sum(axis=1) if i > 0 else None
    ax.bar(isolates, [1] * len(isolates), label=cluster, color=data[cluster], bottom=bottom)

# Personnaliser le graphique
plt.xticks(rotation=90, fontsize=8)
plt.legend(title="Clusters", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.ylabel("Proportion")
plt.title("Clustering Results")
plt.tight_layout()

# Enregistrer le graphique sous forme d'image
output_path = "/media/name/ESTEL3/TREE_CIBIG/clustering_results.png"
plt.savefig(output_path)
print(f"Graphique enregistré dans : {output_path}")

