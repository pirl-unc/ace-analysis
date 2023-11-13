#!/bin/bash

#SBATCH --time=4:00:00
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -n 32 # Number of cores
#SBATCH --mem=128G
#SBATCH --job-name=generate-configurations
#SBATCH -o 03_generate_configurations_slurm.out

export OMP_NUM_THREADS=1
python 03_generate_configurations.py
