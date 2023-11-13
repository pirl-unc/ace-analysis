#!/bin/bash

#SBATCH --time=2:00:00
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -n 64 # Number of cores
#SBATCH --mem=256G
#SBATCH --job-name=generate-configurations
#SBATCH -o 03_generate_configurations_slurm.out

export OMP_NUM_THREADS=1
python 03_generate_configurations.py
