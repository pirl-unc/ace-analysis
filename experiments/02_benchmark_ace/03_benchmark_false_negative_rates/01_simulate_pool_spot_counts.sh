#!/bin/bash

#SBATCH --time=4:00:00
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -n 64 # Number of cores
#SBATCH --mem=256G
#SBATCH --job-name=simulate-pool-spot-counts
#SBATCH -o 01_simulate_pool_spot_counts_slurm.out

export OMP_NUM_THREADS=1
python 01_simulate_pool_spot_counts.py
