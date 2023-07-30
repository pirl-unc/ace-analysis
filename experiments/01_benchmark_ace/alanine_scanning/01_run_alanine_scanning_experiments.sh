#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -n 24 # Number of cores
#SBATCH --mem=96G
#SBATCH --job-name=run-alanine-scanning-experiments
#SBATCH -o run_alanine_scanning_experiments_slurm-%A_%a.out

export OMP_NUM_THREADS=1
python 01_run_alanine_scanning_experiments.py
