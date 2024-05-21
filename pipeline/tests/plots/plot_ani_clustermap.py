from __future__ import annotations

import itertools
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as hc
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from seaborn.matrix import ClusterGrid


def run(
    ani_df: pd.DataFrame,
    outdir: Path,
    dendrogram_ratio: float = 0.15,
    cmap_colors: list[str] | None = None,
    cmap_gamma: float = 1.0,
    cbar_pos: tuple[float, float, float, float] = (0.02, 0.8, 0.05, 0.18),
) -> None:
    # Read ANI matrix
    all_values = itertools.chain.from_iterable(ani_df.values)
    min_ani = min(filter(lambda v: v != 0, all_values))

    # Hierarchical clustering ANI matrix
    linkage = hc.linkage(ani_df.values, method="average")

    # Draw ANI clustermap
    cmap_colors = ["lime", "yellow", "red"] if cmap_colors is None else cmap_colors
    mycmap = LinearSegmentedColormap.from_list(
        "mycmap", colors=cmap_colors, gamma=cmap_gamma
    )
    mycmap.set_under("lightgrey")

    g: ClusterGrid = sns.clustermap(
        data=np.floor(ani_df * 10) / 10,
        # method="average",
        col_linkage=linkage,
        row_linkage=linkage,
        figsize=(max(len(ani_df) / 5, 10), max(len(ani_df) / 5, 10)),
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
    ani_clustermap_file = outdir / "ANIclustermap.png"
    plt.savefig(ani_clustermap_file, dpi=600)
    plt.savefig(ani_clustermap_file.with_suffix(".svg"), dpi=600)
    plt.close()


if __name__ == "__main__":
    ani_df = pd.read_csv(Path("ani_pairwise_values.csv"), index_col='query')
    ani_df.index.name = "Genome"
    ani_df.drop(columns=['max_value', 'max_column'], inplace=True)
    df_large = pd.DataFrame(ani_df.values, columns=list(ani_df.columns), index=list(ani_df.index))
    # ani_df = ani_df.iloc[:2, :2]
    run(df_large, Path.cwd())
