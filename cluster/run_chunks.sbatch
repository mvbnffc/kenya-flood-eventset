#!/bin/bash
# Simulate the event pool in chunks: one array task per chunk, each running 05b_simulate_chunk.R.
# Submit from the repo root with the array size 05a prints:  sbatch --array=1-N%100 cluster/run_chunks.sbatch
#SBATCH --job-name=kenya_flood_05b
#SBATCH --partition=Short
# --array=1-N%100 is given on the sbatch command line (N from 05a), as N changes with the config.
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=cluster/logs/chunk_%a_%A.out
#SBATCH --error=cluster/logs/chunk_%a_%A.err

set -euo pipefail
export OMP_NUM_THREADS=1   # one CPU per task
/usr/bin/Rscript scripts/05b_simulate_chunk.R "$SLURM_ARRAY_TASK_ID"
