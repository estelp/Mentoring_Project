import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from PIL import ImageColor

# Charger le fichier Excel
file_path = '/media/name/ESTEL3/TREE_CIBIG/K_file.xlsx'
data = pd.read_excel(file_path)

# Vérifier les premières lignes des données
print(data.head())

# Supposons que la première colonne contient les identifiants des isolats
# Nous excluons cette colonne pour ne conserver que les couleurs associées aux différentes valeurs de K
color_data = data.iloc[:, 1:]  # Colonnes restantes pour les couleurs

# Fonction pour mélanger les couleurs
def mix_colors(colors):
    # Convertir les noms de couleurs en RGBA (composants compris entre 0 et 1)
    rgba_colors = [ImageColor.getcolor(color, "RGBA") for color in colors]
    
    # Normaliser les valeurs RGB dans la plage [0, 1] pour le calcul de la moyenne
    rgba_colors_normalized = [(r / 255, g / 255, b / 255, a / 255) for r, g, b, a in rgba_colors]
    
    # Calculer la couleur moyenne (moyenne des canaux RGBA)
    avg_color = tuple(sum(c[i] for c in rgba_colors_normalized) / len(rgba_colors_normalized) for i in range(4))
    
    # Retourner la couleur moyenne sous forme de code hexadécimal (en prenant uniquement RGB)
    return mcolors.to_hex(avg_color[:3])

# Calculer les couleurs mélangées pour chaque ligne
mixed_colors = [mix_colors(row) for row in color_data.values]

# Création du diagramme unique
fig, ax = plt.subplots(figsize=(3, 8))

# Ajouter les sections avec les couleurs mélangées
for i, color in enumerate(mixed_colors):
    ax.bar(0, 1, bottom=i, color=color, edgecolor='black', width=0.5)
    ax.text(0.6, i + 0.5, color, va='center', fontsize=8)

# Ajustements visuels
ax.set_xticks([])
ax.set_yticks([])
ax.set_title('Merged Membership Probabilities')
plt.tight_layout()

# Sauvegarde du graphique dans un fichier
output_path = 'merged_membership_probabilities_plot.png'
plt.savefig(output_path)
print(f"Graphique sauvegardé dans : {output_path}")

