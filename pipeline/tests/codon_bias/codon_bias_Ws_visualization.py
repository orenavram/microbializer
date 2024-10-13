import os
import time

import matplotlib.pyplot as plt
from matplotlib import colors, patches as mpatches
import pandas as pd
import numpy as np
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
import umap.umap_ as umap


def create_2D_TSNE_vector(output_file, WData):
    start_time = time.time()
    """
    input: fileList of highly expressed nucleotide fastas

    output: 2D vector (reduced with t-SNE)
    """
    filepath = os.path.join(output_file, 'genomeIndex')
    fileList = os.listdir(filepath)
    dataframe = pd.DataFrame(WData)
    dataframe = dataframe.transpose()

    reducer = TSNE(n_components=2, random_state=42, perplexity=2)  # 2 components used to reduce to 2 dimensions
    values_reduced = reducer.fit_transform(dataframe)

    # Create the graph
    x = values_reduced[:, 0]
    y = values_reduced[:, 1]
    plt.axis('equal')
    plt.scatter(x, y)
    plt.title("Relative Adaptiveness(W)")

    # Add labels to outliers or all points
    x_mean = x.mean()
    y_mean = y.mean()
    x_std = x.std()
    y_std = y.std()

    MAX_NUMBER_LABELED = 20
    if len(fileList) == len(values_reduced):
        if len(fileList) > MAX_NUMBER_LABELED:
            threshold = 2  # Threshold number of standard deviations to be considered an outlier
            for i in range(len(values_reduced)):
                if abs(x[i] - x_mean) > threshold * x_std or abs(y[i] - y_mean) > threshold * y_std:
                    plt.text(x[i], y[i], f'{fileList[i]}', ha='center', va='bottom')
        else:
            for i in range(len(values_reduced)):
                plt.text(x[i], y[i], f'{fileList[i]}', ha='center', va='bottom')
    else:
        print("Error: Length of fileList and values_reduced don't match.")

    # Make sure output file exists
    if not os.path.exists(output_file):
        os.makedirs(output_file)

    filepath = os.path.join(output_file, 'Relative_Adaptiveness_scatter_plot.png')
    # Save plot to output file
    plt.savefig(filepath)
    plt.close()
    print("Time for t-SNE:", time.time() - start_time)


def create_2D_UMAP_vector(output_file, WData):
    start_time = time.time()
    """
    input: fileList of highly expressed nucleotide fastas

    output: 2D vector (reduced with UMAP) with colored clusters and legend
    """
    filepath = os.path.join(output_file, 'genomeIndex')
    fileList = os.listdir(filepath)
    dataframe = pd.DataFrame(WData)
    dataframe = dataframe.transpose()

    # Perform UMAP dimensionality reduction
    reducer = umap.UMAP(n_components=2, random_state=42)
    values_reduced = reducer.fit_transform(dataframe)

    # Perform K-means clustering to identify clusters
    n_clusters = 5  # Number of clusters
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    cluster_labels = kmeans.fit_predict(values_reduced)

    # Normalize cluster values for consistent colormap
    norm = colors.Normalize(vmin=0, vmax=n_clusters - 1)

    # Create the scatter plot with colored clusters
    scatter = plt.scatter(values_reduced[:, 0], values_reduced[:, 1], c=cluster_labels, cmap='viridis', norm=norm)
    plt.title("Relative Adaptiveness(W)")

    # Create custom legend with colors matching the scatter plot
    cluster_legend = []
    for cluster in range(n_clusters):
        cluster_color = scatter.cmap(norm(cluster))
        legend_patch = mpatches.Patch(color=cluster_color, label=f'Cluster {cluster + 1}')
        cluster_legend.append(legend_patch)

    # Add legend to the plot
    plt.legend(handles=cluster_legend)

    # Add labels to outliers or all points
    x_mean = values_reduced[:, 0].mean()
    y_mean = values_reduced[:, 1].mean()
    x_std = values_reduced[:, 0].std()
    y_std = values_reduced[:, 1].std()

    MAX_NUMBER_LABELED = 20
    if len(fileList) == len(values_reduced):
        if len(fileList) > MAX_NUMBER_LABELED:
            threshold = 2  # Threshold number of standard deviations to be considered an outlier
            for i in range(len(values_reduced)):
                if abs(values_reduced[i, 0] - x_mean) > threshold * x_std or abs(
                        values_reduced[i, 1] - y_mean) > threshold * y_std:
                    plt.text(values_reduced[i, 0], values_reduced[i, 1], f'{fileList[i]}', ha='center', va='bottom')
        else:
            for i in range(len(values_reduced)):
                plt.text(values_reduced[i, 0], values_reduced[i, 1], f'{fileList[i]}', ha='center', va='bottom')
    else:
        print("Error: Length of fileList and values_reduced don't match.")

    # Write labels sorted by clusters to a file
    labels_file = os.path.join(output_file, 'labels_by_cluster.txt')
    with open(labels_file, 'w') as f:
        for cluster in range(n_clusters):
            cluster_indices = np.where(cluster_labels == cluster)[0]
            cluster_labels_sorted = [fileList[i] for i in cluster_indices]
            f.write(f'Cluster {cluster + 1}:\n')
            f.write('\n'.join(cluster_labels_sorted))
            f.write('\n\n')

    # Make sure output file exists
    if not os.path.exists(output_file):
        os.makedirs(output_file)

    filepath = os.path.join(output_file, 'Relative_Adaptiveness_UMAP.png')
    # Save plot to output file
    plt.savefig(filepath)
    plt.close()
    print("Time for UMAP:", time.time() - start_time)
