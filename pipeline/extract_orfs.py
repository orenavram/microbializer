from Bio import SeqIO
import subprocess


# module load prodigal/prodigal-2.6.3
def find_genes(genome, output_path, log_path):
    '''
    input:path to fasta file with prokaryotic genome to be analized
    output: protein-coding gene prediction for input genome
    '''
    segments = list(SeqIO.parse(genome, 'fasta'))
    length = sum(len(segment) for segment in segments)
    if length <= 20000:
        mode = 'meta'
    else:
        mode = 'single'
    cmd = f'prodigal -i {genome}  -d {output_path} -p {mode}'
    #cmd = 'prodigal -i {input_file}  -d {output_file} -o {log_file} -p {Specify_mode}'.format(input_file=genome,
    #                                                                                         output_file=output_path,
    #                                                                                         log_file = log_path,
    #                                                                                         Specify_mode=mode)
    logger.info(f'Starting prodigal. Executed command is:\n{cmd}')
    subprocess.run(cmd, shell=True)

    #TODO: move this part to a preprocessing step
    #remove pipes ("|") from file to avoid queue bugs downstream
    with open(output_path) as f:
        content = f.read()
    content = content.replace("|", "_")
    with open(output_path, 'w') as f:
        f.write(content)


if __name__ == '__main__':
    import logging
    logger = logging.getLogger('main')
    from sys import argv
    logger.info(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

    import argparse
    logger = logging.getLogger('main')

    parser = argparse.ArgumentParser()
    parser.add_argument('genome_path', help='path to fasta genome file')
    parser.add_argument('output_path', help='path to output file')
    parser.add_argument('log_path', help='path to translated sequences')
    args = parser.parse_args()

    find_genes(args.genome_path, args.output_path, args.log_path)

