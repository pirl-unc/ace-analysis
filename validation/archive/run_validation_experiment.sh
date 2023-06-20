#!/bin/bash

#SBATCH --time=72:00:00
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -n 64 # Number of cores
#SBATCH --mem=128G
#SBATCH --job-name=run-validation-experiment
#SBATCH -o run_validation_experiment_slurm.out

python run_validation_experiment.py
