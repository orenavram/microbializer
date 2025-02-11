import itertools
import argparse
from pathlib import Path
import sys
import pandas as pd
import re
from ete3 import orthoxml

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent.parent))

from pipeline.auxiliaries.run_step_utils import add_default_step_args, run_step
from pipeline.auxiliaries.general_utils import str_to_bool
from pipeline.auxiliaries import consts

ORPHANS_FILENAME_GENOME_NAME_PATTERN = re.compile('(.+)_orphans.txt')


def parse_species_name(logger, name, qfo_benchmark=False):
    if not qfo_benchmark:
        return name, None, None

    try:
        database_name, ncbi_tax_id, _ = name.split('_')
        return None, ncbi_tax_id, database_name
    except Exception:
        logger.exception(f"Failed parsing species name {name}, using full name")
        return name, None, None


def parse_gene_name(logger, name, qfo_benchmark=False):
    if not qfo_benchmark:
        return name

    try:
        prot_id = name.split('|')[1]
        return prot_id
    except Exception:
        logger.exception(f"Failed parsing gene name {name}, using full name")
        return name


def fix_orthoxml_output_file(orthoxml_file_path):
    with open(orthoxml_file_path, 'r') as oxml_file:
        content = oxml_file.read()

    # remove b'' prefixes
    pattern = re.compile(r'b\'\"([\w.:|-]+)\"\'')
    fixed_content = pattern.sub(r'"\1"', content)

    # remove 'ortho:' prefixes
    fixed_content = fixed_content.replace('ortho:', '')

    # add xmlns attribute to root xml tag
    fixed_content = fixed_content.replace('<orthoXML ', r'<orthoXML xmlns="http://orthoXML.org/2011/" ')

    with open(orthoxml_file_path, 'w') as oxml_file:
        oxml_file.write(fixed_content)


def build_orthoxml_and_tsv_output(logger, all_clusters_df, output_dir, qfo_benchmark=False):
    gene_name_to_gene_id = {}
    gene_name_to_gene_prot_id = {}

    # Creates an empty orthoXML object
    oxml = orthoxml.orthoXML(version="0.3", origin="Microbializer", originVersion="1")

    # Add the species (and their genes) to the orthoXML document
    next_id = 1
    for strain_column in all_clusters_df.columns[1:]:
        species_name, ncbi_tax_id, database_name = parse_species_name(logger, strain_column, qfo_benchmark)
        species_xml = orthoxml.species(name=species_name, NCBITaxId=ncbi_tax_id)
        database_xml = orthoxml.database(name=database_name)
        genes_xml = orthoxml.genes()

        all_strain_genes = all_clusters_df[strain_column]
        for og_index, og_strain_genes in all_strain_genes.items():
            if pd.isna(og_strain_genes):
                continue

            try:
                for gene_name in og_strain_genes.split(';'):
                    prot_id = parse_gene_name(logger, gene_name, qfo_benchmark)
                    gene_xml = orthoxml.gene(id=str(next_id), protId=prot_id)
                    genes_xml.add_gene(gene_xml)

                    gene_name_to_gene_id[gene_name] = next_id
                    gene_name_to_gene_prot_id[gene_name] = prot_id

                    next_id += 1
            except Exception as e:
                raise BrokenPipeError(f"Failed parsing gene names for column {strain_column} and OG {og_index}: {e}")

        database_xml.set_genes(genes_xml)
        species_xml.add_database(database_xml)
        oxml.add_species(species_xml)

    # Add an ortho group container to the orthoXML document
    groups_xml = orthoxml.groups()
    oxml.set_groups(groups_xml)

    # Add ortholog groups to the orthoXML document + construct a list of all ortholog groups
    ortholog_groups = []
    for index, group_row in all_clusters_df.iterrows():
        og_xml = orthoxml.group(id=group_row['OG_name'])
        og_list = []
        for strain_genes in group_row[1:]:
            if pd.isna(strain_genes):
                continue

            strain_genes = strain_genes.split(';')
            if len(strain_genes) == 1:
                gene = strain_genes[0]
                gene_ref_xml = orthoxml.geneRef(gene_name_to_gene_id[gene])
                og_xml.add_geneRef(gene_ref_xml)
                og_list.append([gene_name_to_gene_prot_id[gene]])
            else:
                paralog_group_xml = orthoxml.group()
                paralog_group_list = []
                for gene in strain_genes:
                    gene_ref_xml = orthoxml.geneRef(gene_name_to_gene_id[gene])
                    paralog_group_xml.add_geneRef(gene_ref_xml)
                    paralog_group_list.append(gene_name_to_gene_prot_id[gene])
                og_xml.add_paralogGroup(paralog_group_xml)
                og_list.append(paralog_group_list)

        groups_xml.add_orthologGroup(og_xml)

        if consts.OUTPUT_TSV_OF_ORTHOLOGS_PAIRS:
            ortholog_groups.append(og_list)

    # export orthoXML document to output_file
    orthoxml_output_file_path = output_dir / 'orthogroups.orthoxml'
    with open(orthoxml_output_file_path, "w") as oxml_file:
        oxml.export(oxml_file, level=0)

    fix_orthoxml_output_file(orthoxml_output_file_path)
    logger.info(f'Created orthoxml output at {orthoxml_output_file_path}')

    # write to tsv output all ortholog pairs
    if consts.OUTPUT_TSV_OF_ORTHOLOGS_PAIRS:
        tsv_output_file_path = output_dir / 'ortholog_pairs.tsv'
        with open(tsv_output_file_path, "w") as tsv_file:
            for og_list in ortholog_groups:
                for strain1_genes, strain2_genes in itertools.combinations(og_list, 2):
                    for strain1_gene, strain2_gene in itertools.product(strain1_genes, strain2_genes):
                        tsv_file.write(f'{strain1_gene}\t{strain2_gene}\n')
        logger.info(f'Created tsv orthologs pairs output at {tsv_output_file_path}')


def create_phyletic_pattern(logger, orthogroups_df, output_dir):
    phyletic_patterns_str = ''
    strain_names = list(orthogroups_df.columns[1:])
    for strain_name in strain_names:
        phyletic_pattern = ''.join(pd.notnull(orthogroups_df[strain_name]).astype(int).astype(str))
        phyletic_patterns_str += f'>{strain_name}\n{phyletic_pattern}\n'

    phyletic_patterns_path = output_dir / 'phyletic_pattern.fas'
    with open(phyletic_patterns_path, 'w') as f:
        f.write(phyletic_patterns_str)
    logger.info(f'Created phyletic pattern at {phyletic_patterns_path}')


def create_orthogroups_variations(logger, orthologs_table_path, output_dir, qfo_benchmark):
    orthogroups_df = pd.read_csv(orthologs_table_path)
    create_phyletic_pattern(logger, orthogroups_df, output_dir)
    build_orthoxml_and_tsv_output(logger, orthogroups_df, output_dir, qfo_benchmark)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('orthologs_table_path', type=Path, help='path to the orthologs table (input)')
    parser.add_argument('output_dir', type=Path, help='path to the output dir')
    parser.add_argument('--qfo_benchmark', help='whether the output OrthoXml should be in QfO benchmark format',
                        type=str_to_bool)
    add_default_step_args(parser)
    args = parser.parse_args()

    run_step(args, create_orthogroups_variations, args.orthologs_table_path, args.output_dir, args.qfo_benchmark)
