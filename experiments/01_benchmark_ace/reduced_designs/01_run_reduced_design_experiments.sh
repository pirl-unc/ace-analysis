#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -n 64 # Number of cores
#SBATCH --mem=256G
#SBATCH --job-name=run-reduced-design-experiments
#SBATCH -o 01_run_reduced_design_experiments_slurm-%A_%a.out

export OMP_NUM_THREADS=1
python 01_run_reduced_design_experiments.py
