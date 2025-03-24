import os
import shutil

EXTRACTED_DIR = r"C:\repos\microbializer\benchmark_data\extracted"
ORFS_DIR = r"C:\repos\microbializer\benchmark_data\orfs"

os.makedirs(ORFS_DIR, exist_ok=True)
for folder in os.listdir(EXTRACTED_DIR):
    # Construct full folder path
    folder_path = os.path.join(EXTRACTED_DIR, folder)
    orfs_file_path = os.path.join(folder_path, 'ncbi_dataset', 'data', folder, 'cds_from_genomic.fna')

    # Copy orfs file to ORFS_DIR
    shutil.copy(orfs_file_path, os.path.join(ORFS_DIR, f'{folder}.fna'))
