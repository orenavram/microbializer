#!/bin/bash
# Submit the weekly Microbializer report as a SLURM job.
# Usage: ./submit_weekly_report.sh
#
# Schedule weekly job (crontab, every Monday at 08:00. Scheduled in powerslurm-login.tau.ac.il, user microbi):
#   0 8 * * 1 /bin/bash /microbializer/scripts/submit_weekly_report.sh

set -euo pipefail

LOGS_DIR="/microbializer/weekly_report_logs"
SLURM_SCRIPT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/weekly_report.slurm"

mkdir -p "${LOGS_DIR}"

sbatch "${SLURM_SCRIPT}"
