import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.cluster import KMeans, DBSCAN
from sklearn.metrics import silhouette_score

# Check for input arguments
if len(sys.argv) != 3:
    print("Usage: python clustering_analysis.py <PCA_FILE> <OUTPUT_DIR>")
    sys.exit(1)

# Input and output paths
pca_file = sys.argv[1]
output_dir = sys.argv[2]

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Load PCA data
try:
    print(f"Loading PCA data from {pca_file}...")
    data = pd.read_csv(pca_file, sep=r'\s+', header=None)
    data.columns = ['FID', 'IID', 'PC1', 'PC2', 'PC3']
    X = data[['PC1', 'PC2', 'PC3']]
except Exception as e:
    print(f"Error loading PCA data: {e}")
    sys.exit(1)

# Determine optimal number of clusters using the elbow method
print("Calculating optimal number of clusters (Elbow Method)...")
inertia = []
k_values = range(1, 10)
for k in k_values:
    kmeans = KMeans(n_clusters=k, random_state=42, n_init='auto')
    kmeans.fit(X)
    inertia.append(kmeans.inertia_)

plt.figure(figsize=(8, 5))
plt.plot(k_values, inertia, marker='o')
plt.title('Elbow Method')
plt.xlabel('Number of Clusters (k)')
plt.ylabel('Inertia')
elbow_plot_path = os.path.join(output_dir, 'elbow_method.png')
plt.savefig(elbow_plot_path)
print(f"Elbow Method plot saved to {elbow_plot_path}")

# Apply k-means clustering with k=3 (or adjust based on elbow results)
print("Applying k-means clustering...")
kmeans = KMeans(n_clusters=3, random_state=42, n_init='auto')
data['Cluster'] = kmeans.fit_predict(X)

# Sauvegarder les résultats avec les clusters dans un fichier CSV
output_csv = os.path.join(output_dir, 'clustered_data.csv')
data.to_csv(output_csv, index=False)
print(f"Données annotées avec clusters sauvegardées dans {output_csv}")


# Visualize k-means clusters in 3D
print("Generating 3D k-means plot...")
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
scatter = ax.scatter(data['PC1'], data['PC2'], data['PC3'], c=data['Cluster'], cmap='viridis', s=100)
plt.colorbar(scatter, ax=ax)
ax.set_title('k-means Clustering (k=3)')
ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
ax.set_zlabel('PC3')
kmeans_plot_path = os.path.join(output_dir, 'kmeans_pca_plot_3d.png')
plt.savefig(kmeans_plot_path)
print(f"3D k-means plot saved to {kmeans_plot_path}")

# Apply DBSCAN clustering
print("Applying DBSCAN clustering...")
dbscan = DBSCAN(eps=0.1, min_samples=3)
data['DBSCAN_Cluster'] = dbscan.fit_predict(X)

# Visualize DBSCAN clusters in 3D
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
scatter = ax.scatter(data['PC1'], data['PC2'], data['PC3'], c=data['DBSCAN_Cluster'], cmap='plasma', s=100)
plt.colorbar(scatter, ax=ax)
ax.set_title('DBSCAN Clustering')
ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
ax.set_zlabel('PC3')
dbscan_plot_path = os.path.join(output_dir, 'dbscan_pca_plot_3d.png')
plt.savefig(dbscan_plot_path)
print(f"3D DBSCAN plot saved to {dbscan_plot_path}")

# Calculate silhouette scores for different k
print("Calculating silhouette scores...")
silhouette_scores = []
for k in range(2, 10):
    kmeans = KMeans(n_clusters=k, random_state=42, n_init='auto')
    labels = kmeans.fit_predict(X)
    silhouette_scores.append(silhouette_score(X, labels))

plt.figure(figsize=(8, 5))
plt.plot(range(2, 10), silhouette_scores, marker='o')
plt.title('Silhouette Scores')
plt.xlabel('Number of Clusters (k)')
plt.ylabel('Silhouette Score')
silhouette_plot_path = os.path.join(output_dir, 'silhouette_scores.png')
plt.savefig(silhouette_plot_path)
print(f"Silhouette Scores plot saved to {silhouette_plot_path}")

print("Clustering analysis completed.")

