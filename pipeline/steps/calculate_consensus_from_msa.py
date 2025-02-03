from sys import argv
import argparse
import sys
import traceback
from Bio import SeqIO
from pathlib import Path
from collections import Counter

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent))

from auxiliaries.pipeline_auxiliaries import get_job_logger, add_default_step_args


def calc_consensus(logger, og_path, output_dir):
    og_name = og_path.stem
    consensus_faa_path = output_dir / f"{og_name}.faa"

    if consensus_faa_path.exists():
        return

    # Read sequences into a list
    sequences = [str(record.seq) for record in SeqIO.parse(og_path, "fasta")]
    number_of_sequences = len(sequences)
    alignment_length = len(sequences[0])  # All sequences should be the same length

    consensus = []
    for i in range(alignment_length):
        column = [seq[i] for seq in sequences]  # Extract the column (all bases at position i)

        counts = Counter(column)
        consensus_base, consensus_count = counts.most_common(1)[0]

        # Handle gaps (`-`)
        if consensus_base == '-':
            if consensus_count / number_of_sequences > 0.5:
                continue  # Skip column
            else:
                consensus_base = counts.most_common(2)[1][0]  # Find the most common base that is not a gap

        consensus.append(consensus_base)

    with open(consensus_faa_path, "w") as f:
        f.write(f">{og_name}_consensus\n{''.join(consensus)}\n")

    # og_tmp_dir = os.path.join(output_dir, f'{og_name}_tmp')
    # os.makedirs(og_tmp_dir, exist_ok=True)
    #
    # msa_stockholm_path = os.path.join(og_tmp_dir, f"{og_name}.sto")
    # AlignIO.convert(og_path, "fasta", msa_stockholm_path, "stockholm")
    #
    # msa_db_path = os.path.join(og_tmp_dir, f"{og_name}_msaDb")
    # profile_db_path = os.path.join(og_tmp_dir, f"{og_name}_profileDB")
    # consensus_db_path = os.path.join(og_tmp_dir, f"{og_name}_consensusDb")
    # consensus_faa_raw_path = os.path.join(og_tmp_dir, f"{og_name}_consensus.faa")
    #
    # cmds = [f"mmseqs convertmsa {msa_stockholm_path} {msa_db_path} -v 1",
    #         f"mmseqs msa2profile {msa_db_path} {profile_db_path} --match-mode 1 -v 1",
    #         f"mmseqs profile2consensus {profile_db_path} {consensus_db_path} -v 1",
    #         f"mmseqs result2flat {consensus_db_path} {consensus_db_path} {consensus_db_path} {consensus_faa_raw_path} -v 1"]
    #
    # for cmd in cmds:
    #     logger.info(f'Calling: {cmd}')
    #     subprocess.run(cmd, shell=True, check=True)
    #
    # record = SeqIO.parse(consensus_faa_raw_path, 'fasta').__next__()
    # record.id = f"{og_name}_consensus"
    #
    # SeqIO.write(record, consensus_faa_path, 'fasta')
    # shutil.rmtree(og_tmp_dir, ignore_errors=True)

    logger.info(f'Consensus calculation finished. Output written to {consensus_faa_path}')


def calc_consensus_from_all_ogs(logger, job_input_path, output_dir):
    with open(job_input_path, 'r') as f:
        ogs_paths = [line.strip() for line in f]

    for og_path in ogs_paths:
        calc_consensus(logger, Path(og_path), output_dir)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('job_input_path', type=Path, help='path to job input file with OG paths')
    parser.add_argument('output_dir', type=Path, help='path to output dir')
    add_default_step_args(parser)
    args = parser.parse_args()

    logger = get_job_logger(args.logs_dir, args.job_name, args.verbose)

    logger.info(script_run_message)
    try:
        calc_consensus_from_all_ogs(logger, args.job_input_path, args.output_dir)
    except Exception as e:
        logger.exception(f'Error in {Path(__file__).name}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
