import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from matplotlib import colors, patches as mpatches

Ws = r"C:\Users\TalPNB22\OneDrive\Documents\University\Master\Posters\73 ecoli outputs\W_vectors.csv"
OUTPUT_FILE = r'C:\temp\pic.png'

W_vectors = pd.read_csv(Ws, index_col='Genome')

# Perform K-means clustering
genome_count = len(W_vectors)
if genome_count > 75:
    n_clusters = 5
elif genome_count > 30:
    n_clusters = 4
else:
    n_clusters = 3
kmeans = KMeans(n_clusters=n_clusters, random_state=42)
cluster_labels = kmeans.fit_predict(W_vectors.values)

# Perform PCA dimensionality reduction
pca = PCA(n_components=2)
Ws_values_reduced = pca.fit_transform(W_vectors.values)
explained_variance_ratio = pca.explained_variance_ratio_ * 100

# Normalize cluster values for consistent colormap
norm = colors.Normalize(vmin=0, vmax=n_clusters - 1)

# Create the scatter plot with color-coded clusters
x = Ws_values_reduced[:, 0]
y = Ws_values_reduced[:, 1]
plt.axis('equal')
scatter = plt.scatter(x, y, c=cluster_labels, cmap='viridis', alpha=0.4)
plt.title("Relative Adaptiveness (W vectors) of Genomes", fontsize=20, loc='center', wrap=True)
plt.xlabel(f"PC1 ({explained_variance_ratio[0]:.1f}%)", fontsize=15)
plt.ylabel(f"PC2 ({explained_variance_ratio[1]:.1f}%)", fontsize=15)

cluster_legend = []
for cluster in range(n_clusters):
    cluster_color = scatter.cmap(norm(cluster))
    legend_patch = mpatches.Patch(color=cluster_color, label=f'Cluster {cluster + 1}')
    cluster_legend.append(legend_patch)

plt.legend(handles=cluster_legend)

# Save plot to output_dir
plt.tight_layout()
plt.savefig(OUTPUT_FILE, dpi=600)
plt.close()
