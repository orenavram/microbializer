import subprocess
import os

def generateMsa(sequences_file_path, output_file_path):
    """Generating Multiple sequence alignment using MAFFT module  (alignment method = linsi).

    Parameters
    ----------
    sequences_file_path : string
       Path to FASTA format file containing required sequences to align.
    output_file_path: string
        Path to  desired FASTA format file to contain the MSA.

    Returns
    -------
    FASTA file
        original FASTA format file aligned.
    """
    cmd = 'mafft --maxiterate 1000 --localpair {seqs} > {out}'.format(seqs=sequences_file_path, out=output_file_path)
    subprocess.call(cmd, shell=True)

if __name__ == '__main__':
    import logging
    logger = logging.getLogger('main')
    from sys import argv
    if len(argv) < 2:
        logger.error('Usage: msa <sequences_file_path>')
        exit()
    else:
        generateMsa(*argv[1:])


