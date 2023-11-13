#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -n 96 # Number of cores
#SBATCH --mem=384G
#SBATCH --job-name=run-design-exploration
#SBATCH -o 01_run_design_exploration_slurm.out

export OMP_NUM_THREADS=1
python 01_run_design_exploration.py
