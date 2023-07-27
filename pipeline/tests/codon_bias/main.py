import subprocess
import sys
from sys import argv
import argparse
import logging
import os
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import pandas
import numpy as np
import time
import json
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
import umap.umap_ as umap
from matplotlib import colors
from matplotlib import patches as mpatches

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(os.path.dirname(SCRIPT_DIR)))

from auxiliaries.pipeline_auxiliaries import get_job_logger


def getWData(output_file):
    WData = {}
    filepath = output_file + '/genomeIndex/'
    filedir = os.listdir(filepath)
    for filename in filedir: 
        with open(filepath + filename, "r") as file:
            genomeIndex = json.load(file)
            WData[filename] = genomeIndex
    return WData
    

def create_2D_vector(output_file, WData, gene_number):
    """
    input: location of output file
    Dictionary with genome and Condon Adaption Index 
    fileList of highly expressed nucleotide fastas
   
    output: 2D vector (reduced with t-SNE) Clustered with Kmeans
    """
    start_time = time.time()

    # Create output directory if it doesn't exist
    if not os.path.exists(output_file):
        os.makedirs(output_file)

    # Define file paths
    filepath = os.path.join(output_file, 'genomeIndex')
    points_filepath = os.path.join(output_file, 'point_labels.csv')

    # Get list of file names
    fileList = os.listdir(filepath)

    # Perform PCA dimensionality reduction
    dataframe = pandas.DataFrame(WData).transpose()
    values = dataframe.values


    # Perform K-means clustering
    print(gene_number)
    if(gene_number > 75):
        n_clusters = 5
    elif(gene_number > 30):
        n_clusters = 4
    else:
        n_clusters = 3
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    cluster_labels = kmeans.fit_predict(values)

    # Perform PCA
    pca = PCA(n_components=2)
    values_reduced = pca.fit_transform(values)
    
    # Normalize cluster values for consistent colormap
    norm = colors.Normalize(vmin=0, vmax=n_clusters-1)

    # Create the scatter plot with color-coded clusters
    x = values_reduced[:, 0]
    y = values_reduced[:, 1]
    plt.axis('equal')
    scatter = plt.scatter(x, y, c=cluster_labels, cmap='viridis', alpha=0.4)
    plt.title("Relative Adaptiveness(W)")
    
    cluster_legend = []
    for cluster in range(n_clusters):
        cluster_color = scatter.cmap(norm(cluster))
        legend_patch = mpatches.Patch(color=cluster_color, label=f'Cluster {cluster+1}')
        cluster_legend.append(legend_patch)
        
    plt.legend(handles=cluster_legend)


    # Make sure output file exists
    if not os.path.exists(output_file):
        os.makedirs(output_file)

    filepath = os.path.join(output_file, 'Relative_Adaptiveness_scatter_plot.png')
    # Save plot to output file
    plt.savefig(filepath)
    plt.close()

    # Save point labels and coordinates to file
    point_labels_df = pandas.DataFrame({'Genome': fileList, 'X': x, 'Y': y, 'Cluster': cluster_labels+1})
    point_labels_df.to_csv(points_filepath, index=False)

    print("Time for PCA:", time.time() - start_time)


def create_2D_TSNE_vector(output_file, WData):
    start_time = time.time()
    """
    input: fileList of highly expressed nucleotide fastas
   
    output: 2D vector (reduced with t-SNE)
    """
    filepath = os.path.join(output_file, 'genomeIndex')
    fileList = os.listdir(filepath)
    dataframe = pandas.DataFrame(WData)
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
    dataframe = pandas.DataFrame(WData)
    dataframe = dataframe.transpose()


    # Perform UMAP dimensionality reduction
    reducer = umap.UMAP(n_components=2, random_state=42)
    values_reduced = reducer.fit_transform(dataframe)

    # Perform K-means clustering to identify clusters
    n_clusters = 5  # Number of clusters
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    cluster_labels = kmeans.fit_predict(values_reduced)

    # Normalize cluster values for consistent colormap
    norm = colors.Normalize(vmin=0, vmax=n_clusters-1)

    # Create the scatter plot with colored clusters
    scatter = plt.scatter(values_reduced[:, 0], values_reduced[:, 1], c=cluster_labels, cmap='viridis', norm=norm)
    plt.title("Relative Adaptiveness(W)")

    # Create custom legend with colors matching the scatter plot
    cluster_legend = []
    for cluster in range(n_clusters):
        cluster_color = scatter.cmap(norm(cluster))
        legend_patch = mpatches.Patch(color=cluster_color, label=f'Cluster {cluster+1}')
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
                if abs(values_reduced[i, 0] - x_mean) > threshold * x_std or abs(values_reduced[i, 1] - y_mean) > threshold * y_std:
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
            f.write(f'Cluster {cluster+1}:\n')
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


def get_CAI_Data(output_file):
    """
    input:
    location of output file to access individually calculated CAI data
   
    output: Histogram of CAI distribution and table of mean and standard
    deviation for each orthologous group
    """
    
    filepath = output_file + "/OG_CAIs"
    fileList = os.listdir(filepath)
    
    CAI_Data = {}
    
    #Iterate through the directory of OG group CAIs and extract data
    for filename in fileList:
        with open(filepath + "/" + filename, 'r') as file:
            line = file.readline() #Read in the first line
            #Get the mean of the orthologous group
            mean = float(line[line.index("Mean:")+5:line.index("STD:")].strip())
            #Get the standard deviation of the orthologous group
            std = float(line[line.index("STD:")+4:].strip())
            CAI_Data[filename] = {"mean": mean, "std": std}
    
    sorted_CAI_Data = sorted(CAI_Data.items(), key=lambda x: x[1]["mean"], reverse=True)

    entries = [[filename, data["mean"], data["std"]] for filename, data in sorted_CAI_Data]

    header = ["Orthologous Group", "Mean", "Std"]

    entries.insert(0, header)
    data_array = np.array(entries)
    
    np.savetxt(output_file + "/" + "CAI_Table.csv", data_array, delimiter=",", fmt="%s")

    #Create plot parameters
    plt.title("CAI distribution across OGs")
    plt.xlabel('Value')
    plt.ylabel('Frequency')
    plt.axis('auto')
    
    #Isolate the means
    means = [data[1]["mean"] for data in sorted_CAI_Data]
    
    #Plot the means
    plt.hist(means, bins=30)
    
    #Save the Histogram as a png
    filepath = os.path.join(output_file, 'CAI_Histogram.png')
    plt.savefig(filepath)
    plt.close()


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('ORF_dir', help='path to input fasta directory')
    parser.add_argument('OG_dir', help='path to input Orthologous group directory')
    parser.add_argument('output_dir', help='path to output directory')
    parser.add_argument('tmp_dir', help='path to tmp directory')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')

    args = parser.parse_args()

    level = logging.INFO
    logger = get_job_logger(args.logs_dir, level)
    logger.info(script_run_message)


    output_dir = args.output_file
    tmp_dir = "/temp_dir"
    done_files_dir= "/done"
    
    '''
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
    
    #create a done file -- change later
    if not os.path.exists(tmp_dir):
        os.makedirs(done_files_dir)
    
    step_number = '98'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_find_W'
    script_path_one = os.path.join(args.src_dir, 'steps/find_W')
    script_path_two = os.path.join(args.src_dir, 'steps/find_CAI')
    
    orthologs_aa_sequences_dir_path, pipeline_step_tmp_dir = pa.prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt') 
    '''
    
    #1. Calculate Ws
    ORF_fileList = os.listdir(args.ORF_dir) #list of genomes
    #Send a 'batch' to find the W for each 
    for file in ORF_fileList:
        #change logdir
        cmd = f'python /scratch/brown/gdurante/codonBias/get_W.py {args.ORF_dir} {args.OG_dir} {args.HEG_reference_file} {args.output_file} {file} --logs_dir logsdir'
        logger.info('calling batch' + cmd)
        subprocess.run(cmd, shell = True)  
    
    
    #2. Make Graph
    ORF_fileList = os.listdir(args.ORF_dir) #list of genomes
    WData = {}
    WData = getWData(args.output_file)
    create_2D_vector(args.output_file, WData, len(ORF_fileList))
    
    
    
    #3. Calculate CAIS
    NUMBER_OF_JOBS = 10
    NUMBER_PER_JOB = int(len(os.listdir(args.OG_dir)) / (NUMBER_OF_JOBS-1))
    
    for i in range(NUMBER_OF_JOBS):
        start = i *NUMBER_PER_JOB
        stop = (i * NUMBER_PER_JOB) + NUMBER_PER_JOB

        cmd = f'python /scratch/brown/gdurante/codonBias/calculate_cai.py {args.OG_dir} {args.output_file} {start} {stop} --logs_dir {args.logs_dir}'
        logger.info('calling batch' + cmd)
        subprocess.run(cmd, shell = True)
        
    
    #4. Make graph and table     
    get_CAI_Data(args.output_file)
    
    



   
    
    
