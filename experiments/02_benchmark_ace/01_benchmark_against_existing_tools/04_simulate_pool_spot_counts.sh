#!/bin/bash

#SBATCH --time=4:00:00
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -n 32 # Number of cores
#SBATCH --mem=64G
#SBATCH --job-name=simulate-pool-spot-counts
#SBATCH -o 04_simulate_pool_spot_counts_slurm.out

export OMP_NUM_THREADS=1
python 04_simulate_pool_spot_counts.py
