from sys import argv
import argparse
import logging
import os
import sys
import pandas as pd
import re
from ete3 import orthoxml

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger
from auxiliaries import consts


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
        prot_id = name.split('_')[4]
        return prot_id
    except Exception:
        logger.exception(f"Failed parsing gene name {name}, using full name")
        return name


def fix_orthoxml_output_file(orthoxml_file_path):
    with open(orthoxml_file_path, 'r') as oxml_file:
        content = oxml_file.read()

    # remove b'' prefixes
    pattern = re.compile(r'b\'\"([\w-]+)\"\'')
    fixed_content = pattern.sub(r'"\1"', content)

    # remove 'ortho:' prefixes
    fixed_content = fixed_content.replace('ortho:', '')

    # add xmlns attribute to root xml tag
    fixed_content = fixed_content.replace('<orthoXML ', r'<orthoXML xmlns="http://orthoXML.org/2011/" ')

    with open(orthoxml_file_path, 'w') as oxml_file:
        oxml_file.write(fixed_content)


def build_orthoxml_output(logger, og_table_path, output_file, qfo_benchmark=False):
    logger.info(f'Start to build orthoxml output based on {og_table_path}')

    df = pd.read_csv(og_table_path)
    gene_name_to_gene_id = {}

    # Creates an empty orthoXML object
    oxml = orthoxml.orthoXML(version="0.3", origin="Microbializer", originVersion="1")

    # Add the species (and their genes) to the orthoXML document
    next_id = 1
    for strain_column in df.columns[1:]:
        species_name, ncbi_tax_id, database_name = parse_species_name(logger, strain_column, qfo_benchmark)
        species_xml = orthoxml.species(name=species_name, NCBITaxId=ncbi_tax_id)
        database_xml = orthoxml.database(name=database_name)
        genes_xml = orthoxml.genes()

        strain_genes = df[strain_column]
        for gene_name in strain_genes:
            if pd.isna(gene_name):
                continue
            gene_name_to_gene_id[gene_name] = next_id

            prot_id = parse_gene_name(logger, gene_name, qfo_benchmark)
            gene_xml = orthoxml.gene(id=str(next_id), protId=prot_id)
            genes_xml.add_gene(gene_xml)

            next_id += 1

        database_xml.set_genes(genes_xml)
        species_xml.add_database(database_xml)
        oxml.add_species(species_xml)

    # Add an ortho group container to the orthoXML document
    groups_xml = orthoxml.groups()
    oxml.set_groups(groups_xml)

    # Add ortholog groups to the orthoXML document
    for index, group_row in df.iterrows():
        group_xml = orthoxml.group(id=group_row['OG_name'])
        for gene_name in group_row[1:]:
            if pd.isna(gene_name):
                continue
            gene_ref_xml = orthoxml.geneRef(gene_name_to_gene_id[gene_name])
            group_xml.add_geneRef(gene_ref_xml)
        groups_xml.add_orthologGroup(group_xml)

    # export orthoXML document to output_file
    with open(output_file, "w") as oxml_file:
        oxml.export(oxml_file, level=0)

    fix_orthoxml_output_file(output_file)


def get_verified_clusters_set(verified_clusters_path):
    return set([os.path.splitext(file)[0] for file in os.listdir(verified_clusters_path)])


def finalize_table(logger, putative_orthologs_path, verified_clusters_path, finalized_table_path, phyletic_patterns_path,
                   orthoxml_path, qfo_benchmark, delimiter):
    verified_clusters_set = get_verified_clusters_set(verified_clusters_path)
    logger.info(f'verified_clusters_set:\n{verified_clusters_set}')
    with open(putative_orthologs_path) as f:
        # OG_name,GCF_000008105,GCF_000006945,GCF_000195995,GCF_000007545,GCF_000009505
        final_orthologs_table_header = f.readline().rstrip()
        # remove "OG_name," from header
        # index_of_first_delimiter = putative_orthologs_table_header.index(delimiter)
        # final_orthologs_table_header = putative_orthologs_table_header[index_of_first_delimiter + 1:]
        finalized_table_str = final_orthologs_table_header + '\n'
        strain_names = final_orthologs_table_header.split(delimiter)[1:]  # remove "OG_name," from header
        strain_names_to_phyletic_pattern = dict.fromkeys(strain_names, '')
        og_number = 0
        for line in f:
            first_delimiter_index = line.index(delimiter)
            OG_name = line[:first_delimiter_index]
            group = line.rstrip()[first_delimiter_index + 1:]  # remove temporary name (one of the group members)
            splitted_group = group.split(delimiter)
            if OG_name in verified_clusters_set:
                logger.debug(f'Adding {OG_name} to the final table')
                verified_clusters_set.discard(OG_name)
                finalized_table_str += f'og_{og_number}{delimiter}{group}\n'
                og_number += 1
                # extract phyletic pattern of the current OG
                for i in range(len(splitted_group)):
                    strain_name = strain_names[i]
                    pattern = '0' if splitted_group[i] == '' else '1'
                    strain_names_to_phyletic_pattern[strain_name] += pattern

    with open(finalized_table_path, 'w') as f:
        f.write(finalized_table_str)

    phyletic_patterns_str = ''
    for strain_name in strain_names:
        phyletic_patterns_str += f'>{strain_name}\n{strain_names_to_phyletic_pattern[strain_name]}\n'

    with open(phyletic_patterns_path, 'w') as f:
        f.write(phyletic_patterns_str)

    build_orthoxml_output(logger, finalized_table_path, orthoxml_path, qfo_benchmark)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('putative_orthologs_path', help='path to a file with the putative orthologs sets')
    parser.add_argument('verified_clusters_path', help='path to a directory with the verified clusters')
    parser.add_argument('finalized_table_path', help='path to an output file in which the final table will be written')
    parser.add_argument('phyletic_patterns_path',
                        help='path to an output file in which the phyletic patterns fasta will be written')
    parser.add_argument('orthoxml_path', help='path to an output file in which the OrthoXml will be written')
    parser.add_argument('--delimiter', help='delimiter for the putative orthologs table', default=consts.CSV_DELIMITER)
    parser.add_argument('--qfo_benchmark', help='whether the output OrthoXml should be in QfO benchmark format', action='store_true')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        finalize_table(logger, args.putative_orthologs_path, args.verified_clusters_path,
                       args.finalized_table_path, args.phyletic_patterns_path, args.orthoxml_path,
                       args.qfo_benchmark, args.delimiter)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
