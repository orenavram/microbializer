# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project is

M1CR0B1AL1Z3R is a bioinformatics pipeline for large-scale bacterial genomic analysis. It takes FASTA genome files as input and produces orthogroups, phylogenetic trees, functional annotations, and visualizations. It runs locally (standalone) or on SLURM HPC clusters.

## Running the pipeline

```bash
# Activate the conda environment first
conda activate microbializer

# Simple local run
python pipeline/main.py --input <path_to_genome_fasta_files> --cpus <number_of_cpus>

# Via JSON parameter file
python pipeline/main.py --args_json_path <path_to_json_file>
```

Key parameters: `--identity_cutoff`, `--coverage_cutoff`, `--core_minimal_percentage`, `--bootstrap`, `--outgroup`, `--inputs_fasta_type genomes|orfs`, `--filter_out_plasmids`.

## Setup

```bash
conda env create -f microbializer.yaml
conda activate microbializer

# Download the KEGG database (required for functional annotation)
cd pipeline/data/kegg
gdown 1JFAKeDCjdsLbkhDOPB_pz_mNzOtqV6rP
tar -xzvf prokaryote_database.tar.gz
```

## Architecture

### Entry point and orchestration

`pipeline/main.py` is the sole orchestrator (~1000+ lines). It defines one function per pipeline step (e.g., `step_1_fix_input_files`, `step_2_search_orfs`) and calls them sequentially. Each step checks for a `.done` file before running — this is the checkpointing mechanism that allows interrupted runs to resume.

### Core directories

- `pipeline/steps/` — individual analysis scripts (one per task: ORF prediction, homology search, clustering, phylogeny, etc.). Each is a standalone Python script invoked via subprocess.
- `pipeline/auxiliaries/` — shared utilities:
  - `consts.py` — **environment configuration** (see below)
  - `configuration.py` — `Config` dataclass with all pipeline parameters
  - `run_step_utils.py` — job submission (`submit_job`, `submit_batch`) and `wait_for_results`
  - `general_utils.py`, `logic_utils.py`, `main_utils.py`, `infer_orthogroups_logic.py`
- `pipeline/flask/` — web server integration (shared constants, state machine for the web UI)
- `tests/` — exploratory/benchmark scripts; not a unit test suite

### Environment configuration

**`pipeline/auxiliaries/consts.py`** controls which environment the pipeline targets. Change the `ENV` variable at the top of this file:

```python
ENV = 'standalone'   # local machine, no SLURM
# ENV = 'yair_test'  # SLURM cluster (test)
# ENV = 'yair_prod'  # SLURM cluster (prod)
# ENV = 'windows'    # Windows local
```

`standalone` disables job submission (`USE_JOB_MANAGER = False`); all other envs use SLURM arrays.

### Job execution model

Steps are parallelized by splitting input into batches and submitting SLURM job arrays (or running locally via `ProcessPoolExecutor` in standalone mode). `wait_for_results` polls for `.done` files in the step's tmp directory to track completion. The pipeline writes an `error_file_path` on failure that propagates errors to the caller.

### Key external tools

Prodigal (ORF prediction), MMseqs2 (homology search), MCL (clustering), MAFFT (alignment), HMMER (KEGG annotation), FastANI (ANI computation), RAxML/FastTree/IQtree (phylogeny).

## Tests

The `tests/` directory contains standalone exploratory and benchmark scripts, not a pytest suite. To run individual scripts: `python tests/<subdir>/<script>.py`. There is no global test runner command.
