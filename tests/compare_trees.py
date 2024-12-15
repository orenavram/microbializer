from ete3 import Tree

raxml_path = r"/groups/pupko/yairshimony/test/phylogeny_softwares/dataset_300_salmonela/raxml/final_species_tree.txt"
iqtree_path = r"/groups/pupko/yairshimony/test/phylogeny_softwares/dataset_300_salmonela/iqtree/final_species_tree.txt"

with open(raxml_path) as fp:
    raxml_newick = fp.read().strip()

with open(iqtree_path) as fp:
    iqtree_newick = fp.read().strip()

raxml_tree = Tree(raxml_newick)
iqtree_tee = Tree(iqtree_newick)

rf, rf_max, common_leaves, edges_t1, edges_t2, discarded_edges_t1, discarded_edges_t2 = raxml_tree.robinson_foulds(iqtree_tee, unrooted_trees=True)

pass