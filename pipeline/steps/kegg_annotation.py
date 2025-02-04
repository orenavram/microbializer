import subprocess
import sys
from sys import argv
import argparse
from pathlib import Path
from Bio import SeqIO
import pandas as pd
import traceback

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent))

from auxiliaries.pipeline_auxiliaries import get_job_logger, add_default_step_args, str_to_bool
from auxiliaries import consts


def create_fasta_of_unified_ogs_sequences(logger, og_aa_dir, output_fasta, optimize):
    if optimize:
        records = []
        for og_path in og_aa_dir.iterdir():
            # read only the first sequence from each og
            first_record = SeqIO.parse(og_path, 'fasta').__next__()
            records.append(first_record)

        SeqIO.write(records, output_fasta, 'fasta')
        logger.info(f'Wrote {len(records)} records to {output_fasta} (only 1 gene from each og)')
    else:
        subprocess.run(f'cat {og_aa_dir}/* >> {output_fasta}', shell=True, check=True)
        logger.info(f'Wrote all records to {output_fasta}')


def filter_hmmsearh_output(hmmsearch_output):
    # parse hmmsearch output
    hmmsearch_output_df = pd.read_csv(hmmsearch_output, sep=r'\s+', comment='#', names=consts.HMMSEARCH_OUTPUT_HEADER)
    hmmsearch_output_df = hmmsearch_output_df[['target_name', 'query_name', 'full_score', 'domain_score']]
    hmmsearch_output_df.rename(columns={'target_name': 'gene', 'query_name': 'knum'}, inplace=True)

    # read ko list
    ko_list_df = pd.read_csv(consts.KEGG_KO_LIST_PATH)[['knum', 'threshold', 'score_type']]

    # merge hmmsearch output with ko list
    hmmsearch_output_df = hmmsearch_output_df.merge(ko_list_df, left_on='knum', right_on='knum', how='left')

    # filter out rows with score below threshold
    hmmsearch_output_df = hmmsearch_output_df[
        ((hmmsearch_output_df['score_type'] == 'full') &
         (hmmsearch_output_df['full_score'] >= hmmsearch_output_df['threshold'])) |
        ((hmmsearch_output_df['score_type'] == 'domain') &
         (hmmsearch_output_df['domain_score'] >= hmmsearch_output_df['threshold']))]

    return hmmsearch_output_df


def create_gene_to_og_map(og_table_df):
    melted_df = og_table_df.melt(id_vars='OG_name', value_name='genes', var_name='genome')

    # Split the 'genes' column by ';' and explode the list into separate rows
    melted_df['genes'] = melted_df['genes'].str.split(';')
    exploded_df = melted_df.explode('genes').dropna(subset=['genes'])

    # Select and rename the relevant columns
    result_df = exploded_df[['genes', 'OG_name']].rename(columns={'genes': 'gene'})
    return result_df


def add_kegg_annotations_to_og_table(og_table_path, hmmsearch_output_df):
    og_table_df = pd.read_csv(og_table_path)
    gene_to_og_df = create_gene_to_og_map(og_table_df)
    hmmsearch_output_df = hmmsearch_output_df.merge(gene_to_og_df, left_on='gene', right_on='gene', how='left')

    og_to_knums_df = hmmsearch_output_df.groupby('OG_name')['knum'].apply(
        lambda ko_list: ';'.join(ko_list.dropna().unique()))
    og_to_knums_df = og_to_knums_df.reset_index()

    knum_to_description = pd.read_csv(consts.KEGG_KO_LIST_PATH, index_col='knum')['description'].to_dict()

    def map_descriptions(knum_string):
        knums = knum_string.split(';')
        descriptions = [knum_to_description.get(k, '') for k in knums]  # Map each knum to its description
        return ';'.join(filter(None, descriptions))  # Join non-empty descriptions

    # Apply the function to create a new 'description' column
    og_to_knums_df['knum_description'] = og_to_knums_df['knum'].apply(map_descriptions)

    # Merge the new columns into the original table
    og_table_with_kegg_df = og_table_df[['OG_name']].merge(og_to_knums_df, left_on='OG_name', right_on='OG_name',
                                                           how='left')
    return og_table_with_kegg_df


def kegg_annotation(logger, og_aa_dir, og_table_path, output_dir, output_og_table_path, cpus, optimize):
    if output_og_table_path.exists():
        logger.info(f'{output_og_table_path} already exists. Exiting...')
        return

    unified_ogs_sequences = output_dir / 'unified_ogs_sequences.faa'
    create_fasta_of_unified_ogs_sequences(logger, og_aa_dir, unified_ogs_sequences, optimize)

    # run hmmsearch
    hmmsearch_output = output_dir / 'hmmsearch_output.txt'
    cmd = f'hmmsearch --noali -o /dev/null --cpu {cpus} --tblout {hmmsearch_output} {consts.KEGG_DATABASE_PATH} {unified_ogs_sequences}'
    logger.info(f'Running: {cmd}')
    subprocess.run(cmd, shell=True, check=True)
    logger.info(f'Finished running hmmsearch. Output written to {hmmsearch_output}')

    hmmsearch_output_df = filter_hmmsearh_output(hmmsearch_output)
    filtered_hmmsearch_output = output_dir / 'filtered_hmmsearch_output.csv'
    hmmsearch_output_df.to_csv(filtered_hmmsearch_output, index=False)
    logger.info(f'Wrote filtered hmmsearch output to {filtered_hmmsearch_output}')

    og_table_with_kegg_df = add_kegg_annotations_to_og_table(og_table_path, hmmsearch_output_df)
    og_table_with_kegg_df.to_csv(output_og_table_path, index=False)
    logger.info(f'Wrote og table with kegg annotations to {output_og_table_path}')


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('og_aa_dir', type=Path, help='path to a dir of the amino acid sequences of all ogs')
    parser.add_argument('og_table_path', type=Path, help='path to the og table')
    parser.add_argument('output_dir', type=Path, help='path to the output dir')
    parser.add_argument('output_og_table_path', type=Path, help='path to the output og table with kegg annotations')
    parser.add_argument('cpus', help='number of cpus to use')
    parser.add_argument('--optimize', help='whether to use only 1 gene from each og or all genes', type=str_to_bool)
    add_default_step_args(parser)
    args = parser.parse_args()

    logger = get_job_logger(args.logs_dir, args.job_name, args.verbose)

    logger.info(script_run_message)
    try:
        kegg_annotation(logger, args.og_aa_dir, args.og_table_path, args.output_dir, args.output_og_table_path,
                        args.cpus, args.optimize)
    except Exception as e:
        logger.exception(f'Error in {Path(__file__).name}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
