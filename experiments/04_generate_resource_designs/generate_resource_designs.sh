#!/bin/bash

#SBATCH --time=8:00:00
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -n 32 # Number of cores
#SBATCH --mem=128G
#SBATCH --job-name=generate-resource-designs
#SBATCH -o generate_resource_designs_slurm.out

python generate_resource_designs.py
