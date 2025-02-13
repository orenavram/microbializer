import pandas as pd
from pathlib import Path
from ete3 import Tree


SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR / 'salmonella_300_c_001'

ORIGINAL_ORTHOGROUPS_INFERENCE_ALL_CORE_OGS = DATA_DIR / 'original' / 'tree_all_core_ogs' / 'final_species_tree.newick'
ORIGINAL_ORTHOGROUPS_INFERENCE_1000_CORE_OGS = DATA_DIR / 'original' / 'tree_1000_core_ogs' / 'final_species_tree.newick'
OPTIMIZED_ORTHOGROUPS_INFERENCE_ALL_CORE_OGS = DATA_DIR / 'pseudo_batchsize_mincomparisons' / 'tree_all_core_ogs' / 'final_species_tree.newick'
OPTIMIZED_ORTHOGROUPS_INFERENCE_1000_CORE_OGS = DATA_DIR / 'pseudo_batchsize_mincomparisons' / 'tree_1000_core_ogs' / 'final_species_tree.newick'


def compare_trees(t1, t2, t1_name, t2_name):
    rf, rf_max, leaves_names, splits_t1, splits_t2, discarded_splits_t1, discarded_splits_t2 = t1.robinson_foulds(t2, unrooted_trees=True)

    result = {
        't1_name': t1_name,
        't2_name': t2_name,
        'rf': rf,
        'rf_max': rf_max,
    }

    result_df = pd.DataFrame(result, index=[0])
    return result_df


def main():
    with open(ORIGINAL_ORTHOGROUPS_INFERENCE_ALL_CORE_OGS, "r") as f:
        original_orthogorups_inference_all_core_ogs_tree = Tree(f.readline().strip())

    with open(ORIGINAL_ORTHOGROUPS_INFERENCE_1000_CORE_OGS, "r") as f:
        original_orthogorups_inference_1000_core_ogs_tree = Tree(f.readline().strip())

    with open(OPTIMIZED_ORTHOGROUPS_INFERENCE_ALL_CORE_OGS, "r") as f:
        optimized_orthogorups_inference_all_core_ogs_tree = Tree(f.readline().strip())

    with open(OPTIMIZED_ORTHOGROUPS_INFERENCE_1000_CORE_OGS, "r") as f:
        optimized_orthogorups_inference_1000_core_ogs_tree = Tree(f.readline().strip())

    result1_df = compare_trees(original_orthogorups_inference_all_core_ogs_tree,
                               original_orthogorups_inference_1000_core_ogs_tree,
                               'original_all_core_ogs', 'original_1000_core_ogs')

    result2_df = compare_trees(optimized_orthogorups_inference_all_core_ogs_tree,
                               optimized_orthogorups_inference_1000_core_ogs_tree,
                               'optimized_all_core_ogs', 'optimized_1000_core_ogs')

    result3_df = compare_trees(original_orthogorups_inference_all_core_ogs_tree,
                               optimized_orthogorups_inference_all_core_ogs_tree,
                               'original_all_core_ogs', 'optimized_all_core_ogs')

    result4_df = compare_trees(original_orthogorups_inference_all_core_ogs_tree,
                               optimized_orthogorups_inference_1000_core_ogs_tree,
                               'original_all_core_ogs', 'optimized_1000_core_ogs')

    results_df = pd.concat([result1_df, result2_df, result3_df, result4_df], ignore_index=True)
    results_df.to_csv(DATA_DIR / 'trees_comparison_results.csv', index=False)


if __name__ == "__main__":
    main()
