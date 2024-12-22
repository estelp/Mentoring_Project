import pandas as pd
import matplotlib.pyplot as plt

# Charger le fichier Excel
file_path = '/media/name/ESTEL3/TREE_CIBIG/K_file.xlsx'
data = pd.read_excel(file_path)

# Vérifier les premières lignes des données
print(data.head())

# Supposons que la première colonne contient les identifiants des isolats
# Nous excluons cette colonne pour ne conserver que les couleurs associées aux différentes valeurs de K
isolate_ids = data.iloc[:, 0]  # Première colonne pour les identifiants
color_data = data.iloc[:, 1:]  # Colonnes restantes pour les couleurs

# Création du diagramme
fig, axes = plt.subplots(1, len(color_data.columns), figsize=(15, 40), sharey=True)

# Boucle pour chaque colonne (K)
for i, k in enumerate(color_data.columns):
    colors = color_data[k]  # Les couleurs sont déjà définies dans le fichier
    for j in range(len(colors)):
        axes[i].bar(0, 1, bottom=j, color=colors[j], edgecolor='black', width=0.5)  # Empile les couleurs verticalement
    axes[i].set_title(k)
    axes[i].set_xticks([])
    axes[i].set_yticks([])

# Ajustements visuels
fig.suptitle('Membership Probabilities')
plt.tight_layout()

# Sauvegarde du graphique dans un fichier
output_path = 'membership_probabilities_plot3.png'
plt.savefig(output_path)
print(f"Graphique sauvegardé dans : {output_path}")

