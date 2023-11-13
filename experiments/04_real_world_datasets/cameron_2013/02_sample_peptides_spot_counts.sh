#!/bin/bash

#SBATCH --time=2:00:00
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -n 2 # Number of cores
#SBATCH --mem=16G
#SBATCH --job-name=sample-peptides-spot-counts
#SBATCH -o 02_sample_peptides_spot_counts_slurm.out

export OMP_NUM_THREADS=1
python 02_sample_peptides_spot_counts.py
