import os
import zipfile
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent

DOWNLOADS_DIR = SCRIPT_DIR / 'downloads'
EXTRACTED_DIR = SCRIPT_DIR / 'extracted'
os.makedirs(EXTRACTED_DIR, exist_ok=True)

# Loop through all files in the directory
for file_path in DOWNLOADS_DIR.iterdir():
    if file_path.suffix == '.zip':
        extracted_dir_path = EXTRACTED_DIR / file_path.stem
        if extracted_dir_path.exists():
            print(f"Directory {extracted_dir_path} already exists, skipping extraction.")
            continue
        try:
            # Extract zip file
            with zipfile.ZipFile(file_path, 'r') as zip_ref:
                # Extract to the same directory or specify a different path
                zip_ref.extractall(EXTRACTED_DIR / file_path.stem)
        except Exception:
            print(f"Failed to extract {file_path}")
            continue
