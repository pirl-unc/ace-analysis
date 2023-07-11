#!/bin/bash

#SBATCH --time=72:00:00
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -n 50 # Number of cores
#SBATCH --mem=256G
#SBATCH --job-name=generate-resource-configurations
#SBATCH -o generate_resource_configurations_slurm.out

python generate_resource_configurations.py
