import numpy as np
import dask.dataframe as dd
from pathlib import Path
import subprocess

SCRIPT_DIR = Path(__file__).resolve().parent

RAW_CSV_PATH = r"C:\temp\all_vs_all_raw.m8"
OUTPUT_PARQUET_PATH = r"C:\temp\all_vs_all.filtered"

MMSEQS_OUTPUT_FORMAT = 'query,target,fident,qcov,tcov,evalue,bits'
MMSEQS_OUTPUT_HEADER = MMSEQS_OUTPUT_FORMAT.split(',')
MMSEQS_OUTPUT_COLUMNS_TYPES = {'query': str, 'target': str, 'fident': np.float32, 'qcov': np.float32, 'tcov': np.float32, 'evalue': np.float32, 'bits': np.float32}


def main():

    cmd = f"awk -F'\\t' '$1 != $2' {RAW_CSV_PATH} > {OUTPUT_PARQUET_PATH}"
    subprocess.run(cmd, shell=True)



    return
    m8_df = dd.read_csv(RAW_CSV_PATH, sep='\t', names=MMSEQS_OUTPUT_HEADER,
                        dtype=MMSEQS_OUTPUT_COLUMNS_TYPES)
    m8_df = m8_df[m8_df['query'] != m8_df['target']]
    m8_df['query_genome'] = m8_df['query'].str.split(':').str[0]
    m8_df['target_genome'] = m8_df['target'].str.split(':').str[0]
    m8_df = m8_df[['query', 'query_genome', 'target', 'target_genome', 'bits']]
    m8_df.to_parquet(OUTPUT_PARQUET_PATH)  # Here I always use parquet since it's a huge file


if __name__ == '__main__':
    main()
