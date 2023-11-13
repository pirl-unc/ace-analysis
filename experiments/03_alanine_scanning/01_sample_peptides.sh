#!/bin/bash

#SBATCH --time=4:00:00
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -n 32 # Number of cores
#SBATCH --mem=64G
#SBATCH --job-name=sample-peptides
#SBATCH -o 01_sample_peptides_slurm.out

export OMP_NUM_THREADS=1
python 01_sample_peptides.py
