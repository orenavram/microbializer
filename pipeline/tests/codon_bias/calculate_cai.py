import subprocess
import sys
from sys import argv
import argparse
import logging
import os
from Bio import SeqIO
#from Bio.SeqUtils import CodonUsage
import CodonUsageModified as CodonUsage
#import CAI_Modified as CU
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import pandas
import numpy as np
import time
import pipeline_auxiliaries as pa
import json


def getIndexDict(output_file):
    WDict = {}
    filepath = output_file + '/genomeIndex/'
    filedir = os.listdir(filepath)
    for filename in filedir: 
        with open(filepath + filename, "r") as file:
            index = json.load(file)
            genomeIndex = CodonUsage.CodonAdaptationIndex()
            genomeIndex.set_cai_index(index)
            WDict[filename] = genomeIndex
    return WDict

def calculate_cai(OG_dir, output_file, start, stop, logdir):
    WDict = getIndexDict(output_file)
    fileList = os.listdir(OG_dir)
    
    start = int(start)
    stop = int(stop)

    
    for i in range(stop - start):
        CAI = {}
        CAI_arr = []
        
        if(i + start >= len(fileList)):
            break
        
        with open(OG_dir + "/" + fileList[start + i]) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                CAI[record.id] = (WDict[record.id].cai_for_gene(record.seq))
                CAI_arr.append((WDict[record.id].cai_for_gene(record.seq)))
                
        write_to_file(output_file, CAI, CAI_arr, fileList[start + i])
            
    
def write_to_file(output_file, CAI, CAI_arr, filename):
    #Check to see if file exists
    filepath = (output_file + "/OG_CAIs")
    
    if not os.path.exists(filepath):
        os.makedirs(filepath)
        
    with open(filepath + "/" + filename, 'w') as file:
        line = f"Mean: {np.mean(CAI_arr)} STD: {np.std(CAI_arr)}"
        file.write(line + "\n")

        for value in CAI:
           line = f'{value} {CAI[value]}'
           file.write(line + "\n")
        
        

if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('OG_dir', help='path to input Orthologous group directory')
    parser.add_argument('output_file', help= 'path to output location')
    parser.add_argument('start', help= 'starting index of file (inclusive)')
    parser.add_argument('stop', help= 'stopping index of file')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')


    args = parser.parse_args()

    level = logging.INFO
    logger = pa.get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        calculate_cai(args.OG_dir, args.output_file, args.start, args.stop, args.logs_dir)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
