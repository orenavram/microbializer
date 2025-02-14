from __future__ import annotations

from sys import argv
import itertools
import argparse
import logging
import subprocess
import os
import sys
import traceback
from pathlib import Path
import pandas as pd

import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as hc
from matplotlib.colors import LinearSegmentedColormap
from seaborn.matrix import ClusterGrid

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger


def run_ani(logger, all_genomes_reference_path, output_path,  cpus):
    tmp_results = Path(all_genomes_reference_path).parent
    raw_output_path = os.path.join(tmp_results, 'all_vs_all_raw_output.tsv')

    # No ANI output is reported for a genome pair if ANI value is much below 80% (https://github.com/ParBLiSS/FastANI)
    cmd = f'fastANI --ql {all_genomes_reference_path} --rl {all_genomes_reference_path} -o {raw_output_path} -t {cpus}'
    logger.info(f'Starting fastANI. Executed command is: {cmd}')
    subprocess.run(cmd, shell=True)
    logger.info(f'fastANI finished successfully. Output is saved to {raw_output_path}')

    df = pd.read_csv(raw_output_path, delimiter='\t',
                     names=['query', 'subject', 'ani_value', 'orthologous_segments', 'total_segments'])
    df['query'] = df['query'].apply(lambda path: os.path.splitext(os.path.basename(path))[0])
    df['subject'] = df['subject'].apply(lambda path: os.path.splitext(os.path.basename(path))[0])

    ani_values_df = df.pivot_table(index='query', columns='subject', values='ani_value')
    ani_values_temp_path = os.path.join(tmp_results, 'ani_pairwise_values_temp.csv')
    ani_values_df.to_csv(ani_values_temp_path)
    logger.info(f'ANI temp values were saved to {ani_values_temp_path}')

    plot_ani_clustermap(ani_values_df, Path(output_path))

    if len(ani_values_df) >= 2:
        # Iterate over rows and find max value ignoring diagonal
        max_values = []
        max_values_columns = []
        for i, row in ani_values_df.iterrows():
            row_values = row.drop(i)  # Drop diagonal element
            max_values.append(row_values.max())
            max_values_columns.append(row_values.idxmax())

        ani_values_df['max_value'] = max_values
        ani_values_df['max_column'] = max_values_columns

    ani_values_path = os.path.join(output_path, 'ani_pairwise_values.csv')
    ani_values_df.to_csv(ani_values_path)
    logger.info(f'ANI values were saved to {ani_values_path}')


def plot_ani_clustermap(
    ani_df: pd.DataFrame,
    outdir: Path,
    dendrogram_ratio: float = 0.15,
    cmap_colors: list[str] | None = None,
    cmap_gamma: float = 1.0,
    cbar_pos: tuple[float, float, float, float] = (0.02, 0.8, 0.05, 0.18),
) -> None:
    all_values = itertools.chain.from_iterable(ani_df.values)
    min_ani = min(filter(lambda v: v != 0, all_values))

    cmap_colors = ["lime", "yellow", "red"] if cmap_colors is None else cmap_colors
    mycmap = LinearSegmentedColormap.from_list(
        "mycmap", colors=cmap_colors, gamma=cmap_gamma
    )
    mycmap.set_under("lightgrey")

    figure_size = min(max(len(ani_df) / 5, 10), 60)
    if ani_df.isnull().values.any() or len(ani_df) <= 2:
        # Plot heatmap, since clustermap isn't possible with NaN values
        fig, ax = plt.subplots(figsize=(figure_size, figure_size))
        sns.heatmap(
            data=np.floor(ani_df * 10) / 10,
            annot=len(ani_df) <= 10,
            fmt=".3g",
            cmap=mycmap,
            xticklabels=False,
            yticklabels=True,
            vmin=np.floor(min_ani * 10) / 10,
            vmax=100,
            cbar=True,
            cbar_kws={
                "label": "ANI (%)",
                "orientation": "vertical",
                "spacing": "proportional"
            },
            ax=ax
        )
    else:
        # Hierarchical clustering ANI matrix
        linkage = hc.linkage(ani_df.values, method="average")

        g: ClusterGrid = sns.clustermap(
            data=np.floor(ani_df * 10) / 10,
            # method="average",
            col_linkage=linkage,
            row_linkage=linkage,
            figsize=(figure_size, figure_size),
            annot=len(ani_df) <= 10,
            fmt=".3g",
            cmap=mycmap,
            dendrogram_ratio=dendrogram_ratio,
            xticklabels=False,
            yticklabels=True,
            vmin=np.floor(min_ani * 10) / 10,
            vmax=100,
            cbar=True,
            cbar_pos=cbar_pos,
            cbar_kws={
                "label": "ANI (%)",
                "orientation": "vertical",
                "spacing": "proportional"
            },
            tree_kws={"linewidths": 1.5},
        )

    # Output ANI clustermap figure
    plt.tight_layout()
    ani_map_path = outdir / "ani_map.png"
    plt.savefig(ani_map_path, dpi=600)
    plt.close()
    logger.info(f'ANI clustermap was saved to {ani_map_path}')


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('all_genomes_reference_path', help='path to a file with a list of all genomes paths')
    parser.add_argument('output_path', help='path to the output directory')
    parser.add_argument('--cpus', type=int, default=1, help='Number of CPUs to use')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('--error_file_path', help='path to error file')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        run_ani(logger, args.all_genomes_reference_path, args.output_path, args.cpus)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
