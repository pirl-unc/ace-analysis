#!/bin/bash

#SBATCH --time=4:00:00
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -n 8 # Number of cores
#SBATCH --mem=32G
#SBATCH --job-name=deconvolve-pool-spot-counts
#SBATCH -o 05_deconvolve_pool_spot_counts_slurm.out

export OMP_NUM_THREADS=1
python 05_deconvolve_pool_spot_counts.py
