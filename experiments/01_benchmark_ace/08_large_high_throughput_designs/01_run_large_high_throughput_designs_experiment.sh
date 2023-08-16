#!/bin/bash

#SBATCH --time=96:00:00
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -n 64 # Number of cores
#SBATCH --mem=256G
#SBATCH --job-name=run-large-high-throughput-designs-experiments
#SBATCH -o 01_run_large_high_throughput_designs_experiments_slurm-%A_%a.out

export OMP_NUM_THREADS=1
python 01_run_large_high_throughput_designs_experiment.py
