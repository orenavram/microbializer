import matplotlib
matplotlib.use('Agg')
import os
import regex as re
import subprocess
from Bio import Phylo
import pylab as pb


def extractingCoreGenes(alignments_dir_path, species_amount):
    """Generating a 'species_core_genomes' dictionary - species_core_genomes[specie_a] = core genome of specie_a
        core genomes will contain only genes from alignments containing ALL species (genes aligned = species_amount).

    Parameters
    ----------
    alignments_dir_path : string
       Directory containing multiple FASTA format files of MSA.
       directory should include '/' at the end (ex: '/groups/pupko/rotembar/')
       each MSA file containing orthologs (genes) of different species aligned.
    species_amount : int
       amount of species being analysed in the complete pipeline.

    Returns
    -------
    dictionary
        containing core genome for each specie being analysed in the complete pipeline.
        species_core_genomes[specie_a] = core genome of specie_a
    """
    species_core_genomes = {}
    directory = os.fsencode(alignments_dir_path)
    for file in os.listdir(directory):
        file_name = os.fsdecode(file)
        if file_name.endswith(".fa") or file_name.endswith(".fna"):  # Iterating on FASTA format files only.
            if isCoreGene(alignments_dir_path + file_name, species_amount):   # Adding only core genes to core genomes.
                updateCoreGenomes(alignments_dir_path + file_name, species_core_genomes)
    genomeConcatenating(species_core_genomes)   # Concatenating all aligned genes to one genome per specie.
    return species_core_genomes

def isCoreGene(msa_file_path, species_amount):
    """Returning True if MSA file containing core genes aligned.(ie number pf headers ==  species_amount)

    Parameters
    ----------
    msa_file_path : string
       Path of msa file in FASTA format.
    species_amount : int
       amount of species being analysed in the complete pipeline.

    Returns
    -------
    boolean
        True if the msa containing core genes.
        else, False.
    """
    with open(msa_file_path) as f:
        lines = f.readlines()
    if len(lines) < species_amount * 2:     # minimum check - assuming each gene per line
        return False
    species_counter = 0
    for line in lines:
        if line[0] == '>':  # each header representing specie
            species_counter += 1
    if species_counter == species_amount:
        return True
    else:
        return False

def updateCoreGenomes(msa_file_path, species_core_genomes):
    """Extracting core orthologs (genes) from msa file and updating species core genomes accordingly.

    Parameters
    ----------
    msa_file_path : string
       Path of msa file in FASTA format (sequences aligned == core orthologs (genes)).
    species_core_genomes : dictionary
        Containing list of core genes sequences for each specie being analysed in the complete pipeline.
        species_core_genomes[specie_a] = list of sequences of core genes for specie_a

    Returns
    -------
    None
        updating species_core_genomes (dictionary):
        containing core genome for each specie being analysed in the complete pipeline.
        species_core_genomes[specie_a] = core genome of specie_a.
    """
    with open(msa_file_path) as f:
        lines = f.readlines()
    tmp_seq = []
    new_specie_flag = False
    for line in lines:
        if line[0] == '>':
            if new_specie_flag == True:
                addCoreGene(species_core_genomes, specie_name, ''.join(tmp_seq).rstrip('\n'))  # adding previous specie's gene to genome dict
                new_specie_flag = False
                tmp_seq = []
            specie_name = getSpecieName(line)
            new_specie_flag = True
        else:   # sequence length may be longer than one line
            tmp_seq.append(line)
    addCoreGene(species_core_genomes, specie_name, ''.join(tmp_seq).rstrip('\n'))     # adding last specie's gene to genome dict


def addCoreGene(species_core_genomes, specie_name, gene_seq):
    """Updating species core genomes dictionary:
     adding gene seq for given specie's list of genes

    Parameters
    ----------
    species_core_genomes : dictionary
        Containing list of core genes sequences for each specie being analysed in the complete pipeline.
        species_core_genomes[specie_a] = list of sequences of core genes for specie_a
    specie_name : string
        specie's name (ie dictionary key)
    gene_seq : string
        sequence of core gene

    Returns
    -------
    None
        updating species_core_genomes (dictionary):
        containing list of core genes sequences for each specie being analysed in the complete pipeline.
        species_core_genomes[specie_a] = list of sequences of core genes for specie_a
    """
    if specie_name in species_core_genomes:
        species_core_genomes[specie_name].append(gene_seq)
    else:
        species_core_genomes[specie_name] = [gene_seq]

def getSpecieName(header_line, delimiter='\s+'):
    """Extracting specie's name from header line (based on FASTA format)

    Parameters
    ----------
    header_line : string
    delimiter : string

    Returns
    -------
    string
        specie's name
    """
    res = re.split(delimiter, header_line[1:])
    specie_name = res[0]
    return specie_name

def genomeConcatenating(species_core_genomes):
    """Creating core genomes for each specie by oncatenating each core genes list.

    Parameters
    ----------
       species_core_genomes : dictionary
            Containing list of core genes sequences for each specie being analysed in the complete pipeline.
            species_core_genomes[specie_a] = list of sequences of core genes for specie_a

    Returns
    -------
        None
            updating species_core_genomes (dictionary):
            will contain core genome sequence for each specie being analysed in the complete pipeline.
            species_core_genomes[specie_a] = sequences of core genome for specie_a
        specie's name
    """
    for specie in species_core_genomes:
        species_core_genomes[specie] = ''.join(species_core_genomes[specie])

def generateCoreGenomesMsa(species_core_genomes, core_genomes_msa_file_path='core_msa.fa'):
    """Generating FASTA file containing msa of core genomes for all analysed species in the pipeline.

    Parameters
    ----------
       species_core_genomes : dictionary
            containing core genome for each specie being analysed in the complete pipeline.
            species_core_genomes[specie_a] = core genome of specie_a
        core_genomes_msa_file_path: string
            path for output msa file.

    Returns
    -------
        None
            generating core genomes msa file.
    """

    file = open(core_genomes_msa_file_path, 'w')
    for specie in species_core_genomes.keys():
        file.write('>' + specie)
        file.write('\n')
        file.write(species_core_genomes[specie])
        file.write('\n')
    file.close()


def generatePhylogeneticTree(core_msa_file_path='core_msa.fa', output_file_path='core_genomes_tree'):
    """Generating phylogenetic tree based on core genomes msa using RAXML module.

    Parameters
    ----------
    core_msa_file_path : string
       Path to FASTA format file containing msa of core genomes for all analysed species in the pipeline.
    output_file_path: string
        Path to  desired nexus format file to contain the tree generated.

    Returns
    -------
    None
        generating phylogenetic tree file.
    """
    cmd = 'raxmlHPC -s {msa_path} -n {out} -p {seed} -m {module}'.format(msa_path=core_msa_file_path, out=output_file_path, seed='12345', module='GTRGAMMAI')
    subprocess.call(cmd, shell=True)


def generateTreeFigure(tree_file_directory, output_file_path='tree_figure'):
    """Generating figure based on phylogenetic tree file.

    Parameters
    ----------
    tree_file_directory : string
       Path of file containing newick format tree.
    output_file_path: string
        Path of figure output file.

    Returns
    -------
    None
        generating tree figure.
    """
    tree = Phylo.read(tree_file_directory, 'newick')
    Phylo.draw(tree, do_show=False)
    pb.savefig(output_file_path)


if __name__ == '__main__':
    import logging
    logger = logging.getLogger('main')
    from sys import argv
    if len(argv) < 3:
        logger.error('Usage: phylogenetics_analysis <alignments_dir_path> <species_amount>')
        exit()
    else:
        working_directory = argv[1]
        species_amount = int(argv[2])
        d = extractingCoreGenes(working_directory, species_amount)
        if len(d) == 0:
            logger.error('phylogenetics_analysis - No core genes')
            exit()
        generateCoreGenomesMsa(d, working_directory + 'core_msa.fa')
        generatePhylogeneticTree(working_directory + 'core_msa.fa', 'core_genomes_tree')
        generateTreeFigure(os.getcwd() + '/RAxML_parsimonyTree.core_genomes_tree', working_directory + 'tree_figure')

