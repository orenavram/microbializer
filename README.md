# M1CR0B1AL1Z3R<br>https://microbializer.tau.ac.il

A pipeline for analyzing bacterial genomic sequences. Easily.

<p align="justify">
M1CR0B1AL1Z3R (pronounced: microbializer) was developed to facilitate large-scale microbial genomics analyses. Our aim was to make the analysis easily accessible to the scientific community, addressing the computational challenges of analyzing dozens of bacterial genomes simultaneously. These computations are resource-intensive; therefore, M1CR0B1AL1Z3R leverages high-performance computing clusters to run background processes in parallel. By doing so, it empowers researchers to preprocess and analyze massive genomic datasets without requiring costly computational infrastructure or specialized expertise. The tool generates both visual and textual outputs, ready for publication or further analysis. For more information, visit the <a href="https://microbializer.tau.ac.il/overview">Overview page</a>, and for a running example visit the <a href="https://microbializer.tau.ac.il/gallery">Gallery page</a>.


# Local Run
Although the recommended way to use Microbializer is through the website, we also offer instructions for local installation on your linux machine:
1. Clone the repository
2. Install the microbializer conda environment: 
```angular2html
conda env create -f microbializer.yaml
```
3. Activate the environment: 
```angular2html
conda activate microbializer
```
4. Simple usage: 
```angular2html
python pipeline/main.py --input <path_to_genome_fasta_files> --cpus <number_of_cpus>
```

### Requirements
Microbializer is a computationally intensive pipeline that requires significant computational resources. The recommended way to run it is on an HPC machine or cluster that has many CPU cores and large memory capacity.
The recommended resources are 20 CPUs and 64GB of RAM. However, for small datasets (5 genomes or less), the pipeline can be run with fewer resources, and even with a single CPU core on a standard desktop computer.

### Advanced Usage
For a thorough explanation of the pipeline and the optional parameters, please refer to the website <a href="https://microbializer.tau.ac.il/faq">FAQ page</a>.<br/><br/>
The only required parameter is:
```angular2html
--input <path_to_genome_fasta_files> - Path to a directory/zip/tar.gz containing the input FASTA files.
```
The optional parameters are:
```angular2html
--cpus <number_of_cpus> - Number of CPUs to use. By default uses only 1 core. Recommended to use at least 20 cores for datasets with more than 5 genomes.
--inputs_fasta_type <genomes/orfs> - Type of input FASTA files. Options are: 'genomes' (default) or 'orfs'. In the case of genomes, the first step in the pipeline is ORFs prediction using Prodigal. If the user already has ORFs and would like to use them as the starting point, they should choose ORFs in this parameter, and in that case ORFs prediction is skipped.
--filter_out_plasmids <True/False> - Whether to filter out plasmid contigs/orfs from the input genomes. Default is False. If True, the filtering is done by a simple search for the word "plasmid" in the contig or ORF record name.
--identity_cutoff <int>, --coverage_cutoff <int> - Minimum sequence identity & sequence coverage (protein-level) for homologs detection. These parameters are used in the homology search step, which is the first step in the orthogroup inference step. The homology search is performed using the MMSEQS2 program, which is a fast and sensitive homology search tool. The default values for these parameters are 40(%) sequence identity and 70(%) sequence coverage. These values can be adjusted to increase or decrease the stringency of the homology search.
--core_minimal_percentage <int> - Minimum percent of strains required to consider an orthogroup as part of the core genome. The parameter dictates the inclusion or exclusion of orthogroups in the core proteome (and core genome). By default, this value is set to 100(%) and thus, only orthogroups that contain members of all analyzed genomes are included in the core proteome. However, when bacteria from different orders are analyzed, this strict definition can lead to a very small core proteome. In that case, the core threshold can be lowered. For example, when a 70% threshold is used, orthogroups shared by at least 70% of the analyzed genomes are included in the "core" proteome. In this case, the tree will be inferred using a larger dataset, albeit, with missing values.
--outgroup <name_of_genome> - The outgroup genome is used to root the species tree. By default, no outgroup is used and the produced species tree is unrooted. Alternatively, the user can indicate one of the file names (without the file extension) in the uploaded dataset as an outgroup. In that case, The outgroup genome should be a genome that is phylogenetically distant from the ingroup genomes.
--bootstrap <True/False> - Whether to perform bootstrap analysis for the species tree and display the bootstrap values on it. Default is False.
--add_orphan_genes_to_ogs <True/False> - Default to False, which means that orphan genes (genes that do not belong to any orthogroup) are not included in the orthogroups table (in the result folder 05a_orthogroups). If True, orphan genes are included as orthogroups with single genes. Of note, orthogroups that contain multiple genes of a single genome - are always included in the orthogroups table, regardless of this parameter.
```

An alternative way to supply the parameters is through a json file:
```angular2html
python pipeline/main.py --args_json_path <path_to_json_file>
    
Example of a json file (can include only the parameters you wish to set different than the default):
{
    "input": "<path_to_genome_fasta_files>",
    "cpus": 20,
    "inputs_fasta_type": "genomes",
    "filter_out_plasmids": false,
    "identity_cutoff": 40,
    "coverage_cutoff": 70,
    "core_minimal_percentage": 100,
    "outgroup": null,
    "bootstrap": false,
    "add_orphan_genes_to_ogs": false
}
```

### Outputs
Microbializer provides various outputs divided into sub-folders, as described in the <a href="https://microbializer.tau.ac.il/overview">Overview page</a>.

# Citation 
If you used M1CR0B1AL1Z3R please cite the following papers:
<br/><br/>

M1CR0B1AL1Z3R 2.0: an enhanced web server for comparative analysis of bacterial genomes at scale<br/>
Yair Shimony, Edo Dotan, Elya Wygoda, Naama Wagner, Iris Lyubman, Noa Ecker, Gianna Durante, Gal Mishan, Jeff H Chang, Oren Avram, Tal Pupko<br/> 
Nucleic Acids Research, May 2025; DOI: https://doi.org/10.1093/nar/gkaf413
<br/><br/>

M1CR0B1AL1Z3R - a user-friendly web server for the analysis of large-scale microbial genomics data<br/>
Oren Avram, Dana Rapoport, Shir Portugez, Tal Pupko<br/>
Nucleic Acids Research, May 2019; DOI: https://doi.org/10.1093/nar/gkz423

