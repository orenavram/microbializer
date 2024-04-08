import itertools
import shutil
from sys import argv
import argparse
import logging
import os
import sys
import pandas as pd
import re
from collections import defaultdict
from ete3 import orthoxml

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger
from auxiliaries.logic_auxiliaries import get_strain_name
from auxiliaries import consts

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
        for og_strain_genes in all_strain_genes:
            if pd.isna(og_strain_genes):
                continue

            for gene_name in og_strain_genes.split(';'):
                prot_id = parse_gene_name(logger, gene_name, qfo_benchmark)
                gene_xml = orthoxml.gene(id=str(next_id), protId=prot_id)
                genes_xml.add_gene(gene_xml)

                gene_name_to_gene_id[gene_name] = next_id
                gene_name_to_gene_prot_id[gene_name] = prot_id

                next_id += 1

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
        ortholog_groups.append(og_list)

    # export orthoXML document to output_file
    orthoxml_output_file_path = os.path.join(output_dir, 'orthologs.orthoxml')
    with open(orthoxml_output_file_path, "w") as oxml_file:
        oxml.export(oxml_file, level=0)

    fix_orthoxml_output_file(orthoxml_output_file_path)

    # write to tsv output all ortholog pairs
    tsv_output_file_path = os.path.join(output_dir, 'ortholog_pairs.tsv')
    with open(tsv_output_file_path, "w") as tsv_file:
        for og_list in ortholog_groups:
            for strain1_genes, strain2_genes in itertools.combinations(og_list, 2):
                for strain1_gene, strain2_gene in itertools.product(strain1_genes, strain2_genes):
                    tsv_file.write(f'{strain1_gene}\t{strain2_gene}\n')

    logger.info(f'Created orthoxml and tsv outputs at {orthoxml_output_file_path} and {tsv_output_file_path}')


def finalize_table(logger, orthologs_table_path, finalized_table_path, finalized_table_no_paralogs_path,
                   orphan_genes_dir, qfo_benchmark):
    if orphan_genes_dir:
        orthologs_df = pd.read_csv(orthologs_table_path)

        # Add orphan genes as new OGs
        orphan_clusters = []
        for filename in os.listdir(orphan_genes_dir):
            strain_match_object = ORPHANS_FILENAME_GENOME_NAME_PATTERN.match(filename)
            if not strain_match_object:
                continue
            strain = strain_match_object.group(1)
            with open(os.path.join(orphan_genes_dir, filename)) as orphan_genes_file:
                genes = orphan_genes_file.read().splitlines()
            orphan_clusters.extend([pd.Series({'OG_name': '', strain: gene}) for gene in genes])

        orphan_clusters_df = pd.DataFrame(orphan_clusters)
        all_clusters_df = pd.concat([orthologs_df, orphan_clusters_df], ignore_index=True)
        all_clusters_df['OG_name'] = [f'OG_{i}' for i in range(len(all_clusters_df.index))]
        logger.info(f'Finished adding orphan genes as clusters. OG table now contains {len(all_clusters_df.index)} groups.')

        all_clusters_df.to_csv(finalized_table_path, index=False)

    else:
        shutil.copyfile(orthologs_table_path, finalized_table_path)

    orthologs_df = pd.read_csv(finalized_table_path)
    strain_names = list(orthologs_df.columns[1:])
    output_dir = os.path.dirname(finalized_table_path)

    # Create OG table with only 0 or 1 gene for each genome in each OG
    all_clusters_no_paralogs_df = orthologs_df.copy()
    for index, row in all_clusters_no_paralogs_df.iterrows():
        for strain in strain_names:
            if pd.notnull(row[strain]) and ';' in row[strain]:
                row[strain] = row[strain].split(';')[0]
    all_clusters_no_paralogs_df.to_csv(finalized_table_no_paralogs_path, index=False)

    # Create phyletic pattern
    phyletic_patterns_str = ''
    for strain_name in strain_names:
        phyletic_pattern = ''.join(pd.notnull(orthologs_df[strain_name]).astype(int).astype(str))
        phyletic_patterns_str += f'>{strain_name}\n{phyletic_pattern}\n'

    phyletic_patterns_path = os.path.join(output_dir, 'phyletic_pattern.fas')
    with open(phyletic_patterns_path, 'w') as f:
        f.write(phyletic_patterns_str)
    logger.info(f'Created phyletic pattern at {phyletic_patterns_path}')

    build_orthoxml_and_tsv_output(logger, orthologs_df, output_dir, qfo_benchmark)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('orthologs_table_path', help='path to the orthologs table (input)')
    parser.add_argument('finalized_table_path', help='path to the output orthologs table')
    parser.add_argument('finalized_table_no_paralogs_path', help='path to an output file in which the final table with no paralogs will be written')
    parser.add_argument('--orphan_genes_dir', help='path to a directory with the orphan genes')
    parser.add_argument('--qfo_benchmark', help='whether the output OrthoXml should be in QfO benchmark format', action='store_true')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        finalize_table(logger, args.orthologs_table_path, args.finalized_table_path,
                       args.finalized_table_no_paralogs_path, args.orphan_genes_dir,
                       args.qfo_benchmark)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
