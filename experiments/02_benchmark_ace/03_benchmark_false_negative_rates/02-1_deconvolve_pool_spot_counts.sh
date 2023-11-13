#!/bin/bash

#SBATCH --time=16:00:00
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -n 64 # Number of cores
#SBATCH --mem=256G
#SBATCH --job-name=deconvolve-pool-spot-counts
#SBATCH -o 02-1_deconvolve_pool_spot_counts_slurm.out

export OMP_NUM_THREADS=1
python 02-1_deconvolve_pool_spot_counts.py
