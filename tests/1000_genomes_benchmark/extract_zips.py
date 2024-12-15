import os
import zipfile

DOWNLOADS_DIR = r"C:\repos\microbializer\benchmark_data\downloads"
EXTRACTED_DIR = r"C:\repos\microbializer\benchmark_data\extracted"

os.makedirs(EXTRACTED_DIR, exist_ok=True)
# Loop through all files in the directory
for file in os.listdir(DOWNLOADS_DIR):
    if file.endswith('.zip'):
        # Construct full file path
        file_path = os.path.join(DOWNLOADS_DIR, file)

        # Extract zip file
        with zipfile.ZipFile(file_path, 'r') as zip_ref:
            # Extract to the same directory or specify a different path
            zip_ref.extractall(os.path.join(EXTRACTED_DIR, file[:-4]))
