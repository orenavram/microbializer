import os
import csv
import networkx as nx
import pandas as pd
from pympler import asizeof
import argparse
import logging
import time
from time import sleep
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed, ProcessPoolExecutor

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# HITS_INPUT_DIR = "/home/ai_center/ai_users/yairshimony/microbializer_runs/73_ecoli/outputs_mcl_v2_no_concat_mmseqs_v4_use_parquet/steps_results/05_4_normalize_scores/"
# HITS_INPUT_DIR = "/home/yair/microbializer_runs/4_genomes/outputs_mcl_v2/steps_results/05_4_normalize_scores"
HITS_INPUT_DIR = "/home/ai_center/ai_users/yairshimony/microbializer_runs/salmonella_300/outputs_mcl_v2_mmseqs_v4_use_parquet/steps_results/05_4_normalize_scores/"


def load_file_to_graph(logger, job_input_path):
    """Load edges from CSV files into a graph."""
    start_time = time.time()
    partial_graph = nx.Graph()

    with open(job_input_path, 'r') as f:
        hits_file_paths = f.read().splitlines()

    for hits_path in hits_file_paths:
        with open(hits_path, 'r') as hits_file:
            reader = csv.reader(hits_file)
            next(reader)  # Skip header
            for geneA, geneB, score in reader:
                partial_graph.add_edge(geneA, geneB, weight=round(float(score), 2))

    runtime = time.time() - start_time
    logger.info(f'Partial graph built successfully from file {job_input_path} in {runtime:.2f} seconds')
    return partial_graph


def save_component_to_file(logger, output_dir, component_id, subgraph):
    mcl_file_path = os.path.join(output_dir, f'OG_{component_id}.mcl_input')
    if os.path.exists(mcl_file_path):
        logger.info(f'{mcl_file_path} already exists. Skipping...')
        return

    text_start_time = time.time()
    og_text_for_mcl = ''
    for gene1, gene2, attr in subgraph.edges(data=True):
        og_text_for_mcl += f'{gene1}\t{gene2}\t{attr["weight"]}\n'
    text_runtime = time.time() - text_start_time

    write_start_time = time.time()
    try_index = 1
    while not os.path.exists(mcl_file_path):
        try:
            with open(mcl_file_path, 'w') as mcl_file:
                mcl_file.write(og_text_for_mcl)
        except Exception as e:
            logger.error(f'Error writing {mcl_file_path} (try {try_index}): {e}')
            sleep(1)
            try_index += 1
    write_runtime = time.time() - write_start_time

    logger.info(f'For {mcl_file_path}, text generation took {text_runtime:.2f} seconds, and writing to file took {write_runtime:.2f} seconds.')


def process_files_in_parallel(logger, input_dir, output_dir, num_workers, verbose, parallelize_load_hits, parallelize_write_ogs):
    script_start_time = time.time()

    # Step 1: Build graph from hits files
    logger.info(f"Loading files from {input_dir} and building graph from all files...")
    build_graph_start_time = time.time()
    graph = nx.Graph()
    hits_file_paths = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith('.m8')]
    strains = set([f.split('_vs_')[0] for f in os.listdir(input_dir) if f.endswith('.m8')])
    strains_sorted = sorted(strains)

    if parallelize_load_hits == 'no':
        for hits_file_path in hits_file_paths:
            partial_graph = load_file_to_graph(logger, hits_file_path)
            graph.update(partial_graph)
    else:
        executor_type = ThreadPoolExecutor if parallelize_load_hits == 'threads' else ProcessPoolExecutor

        jobs_inputs_dir = os.path.join(output_dir, 'jobs_inputs')
        os.makedirs(jobs_inputs_dir, exist_ok=True)
        job_index_to_hits_file_paths = defaultdict(list)
        for i, hits_file_path in enumerate(hits_file_paths):
            job_index = i % num_workers
            job_index_to_hits_file_paths[job_index].append(hits_file_path)

        job_input_paths = []
        for job_index, hits_file_paths in job_index_to_hits_file_paths.items():
            job_input_path = os.path.join(jobs_inputs_dir, f'job_{job_index}.txt')
            job_input_paths.append(job_input_path)
            with open(job_input_path, 'w') as job_input_fp:
                job_input_fp.write('\n'.join(hits_file_paths))

        with executor_type(max_workers=num_workers) as executor:
            futures = [executor.submit(load_file_to_graph, logger, job_input_path) for job_input_path in job_input_paths]

            for future in as_completed(futures):
                partial_graph = future.result()
                graph.update(partial_graph)

    build_graph_runtime = time.time() - build_graph_start_time
    logger.info(f"Graph built successfully from all files in {build_graph_runtime:.2f} seconds.")

    if verbose:
        graph_memory = asizeof.asizeof(graph)
        logger.info(f"Total memory size of the graph: {graph_memory} bytes, {graph_memory // 1024} KB, "
                    f"{graph_memory // 1024 ** 2} MB, {graph_memory // 1024 ** 3} GB")

    # Step 2: Find connected components
    components = list(nx.connected_components(graph))
    logger.info(f"Found {len(components)} connected components.")

    # Step 3: Save putative OG table
    logger.info("Writing putative orthologs table to file...")
    write_putative_og_table_start_time = time.time()

    putative_ogs_path = os.path.join(output_dir, 'putative_orthologs_table.csv')
    with open(putative_ogs_path, 'w') as file:
        writer = csv.writer(file)
        writer.writerow(['OG_name'] + strains_sorted)
        for i, component in enumerate(components):
            strain_to_genes = defaultdict(list)
            for gene in component:
                strain = gene.split(':')[0]
                strain_to_genes[strain].append(gene)
            writer.writerow([f'OG_{i}', *[';'.join(sorted(strain_to_genes[strain])) for strain in strains_sorted]])

    write_putative_og_table_runtime = time.time() - write_putative_og_table_start_time
    logger.info(f"Writing putative orthologs table to file {putative_ogs_path} took {write_putative_og_table_runtime:.2f} seconds.")

    # Step 4: Save each putaitve OG to a separate file
    logger.info(f"Writing {len(components)} connected components to files in {output_dir}...")
    write_ogs_start_time = time.time()

    if parallelize_write_ogs == 'no':
        for i, component in enumerate(components):
            subgraph = graph.subgraph(component)
            save_component_to_file(logger, output_dir, i, subgraph)
    else:
        executor_type = ThreadPoolExecutor if parallelize_write_ogs == 'threads' else ProcessPoolExecutor
        with executor_type(max_workers=num_workers) as executor:
            for i, component in enumerate(components):
                subgraph = graph.subgraph(component).copy()
                executor.submit(save_component_to_file, logger, output_dir, i, subgraph)

    write_ogs_runtime = time.time() - write_ogs_start_time
    logger.info(f"Writing connected components to files took {write_ogs_runtime:.2f} seconds.")

    script_runtime = time.time() - script_start_time

    summary = pd.DataFrame({
        'parallelize_load_hits': [parallelize_load_hits],
        'parallelize_write_ogs': [parallelize_write_ogs],
        'build_graph_runtime (minutes)': [build_graph_runtime / 60],
        'write_putative_og_table_runtime (minutes)': [write_putative_og_table_runtime / 60],
        'write_ogs_runtime (minutes)': [write_ogs_runtime / 60],
        'script_runtime (minutes)': [script_runtime / 60],
    })
    summary_path = os.path.join(output_dir, 'summary.csv')
    summary.to_csv(summary_path, index=False)


def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Process MCL input files in parallel.')
    parser.add_argument('--inputs_hits_dir', type=str, required=True, help='Directory containing hits files.')
    parser.add_argument('--num_workers', type=int, default=2, help='Number of worker processes.')
    parser.add_argument('--verbose', action='store_true')
    parser.add_argument('--parallelize_load_hits', type=str, choices=['no', 'threads', 'processes'], default='no')
    parser.add_argument('--parallelize_write_ogs', type=str, choices=['no', 'threads', 'processes'], default='no')
    args = parser.parse_args()

    output_directory = f"{SCRIPT_DIR}/outputs_300_{args.parallelize_load_hits}_{args.parallelize_write_ogs}_{args.num_workers}_load_hits_less_processes/"
    os.makedirs(output_directory, exist_ok=True)

    logger = logging.getLogger(__name__)
    file_handler = logging.FileHandler(os.path.join(output_directory, 'log.txt'), mode='w')
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.setLevel(logging.INFO)

    logger.info(f'Start clustering (parallelize_load_hits={args.parallelize_load_hits}, parallelize_write_ogs={args.parallelize_write_ogs})...')
    process_files_in_parallel(logger, args.inputs_hits_dir, output_directory, args.num_workers, args.verbose,
                              args.parallelize_load_hits, args.parallelize_write_ogs)


if __name__ == "__main__":
    main()
