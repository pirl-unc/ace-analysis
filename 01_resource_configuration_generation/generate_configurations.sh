#!/bin/bash

#SBATCH --time=72:00:00
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -n 50 # Number of cores
#SBATCH --mem=256G
#SBATCH --job-name=generate-configurations
#SBATCH -o generate_configurations_slurm.out

python generate_configurations.py
