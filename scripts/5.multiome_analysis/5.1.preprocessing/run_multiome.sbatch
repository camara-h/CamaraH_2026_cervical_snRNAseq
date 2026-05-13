#!/usr/bin/env bash
#SBATCH --job-name=NK.mlt.ref2020
#SBATCH --partition=medium
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=40
#SBATCH --mem=300G
#SBATCH --output=/home/camarah/storage/Tseng_Neck_Multiome_2026/06_logs/slurm/%x.%j.out
#SBATCH --error=/home/camarah/storage/Tseng_Neck_Multiome_2026/06_logs/slurm/%x.%j.err
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

mkdir -p /home/camarah/storage/Tseng_Neck_Multiome_2026/06_logs/slurm

module purge
module load cellranger/9.0.1
module load cellranger-arc/2.1.0

# conda in batch safe
if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
  source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
  source "$HOME/anaconda3/etc/profile.d/conda.sh"
else
  echo "ERROR: conda.sh not found. Update path." >&2
  exit 1
fi

conda activate snakemake

cd /home/camarah/storage/Tseng_Neck_Multiome_2026/05_scripts

snakemake \
  --cores "${SLURM_CPUS_PER_TASK}" \
  --rerun-incomplete \
  --keep-going \
  --latency-wait 120 \
  -p
