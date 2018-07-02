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
        mode = 'normal'
    cmd = 'prodigal -i {input_file} -d {output_file} -o {log_file} -p {Specify_mode}'.format(input_file=genome,
                                                                                             output_file=output_path,
                                                                                             log_file=log_path,
                                                                                             Specify_mode=mode)
    subprocess.call(cmd.split(), shell=True)
