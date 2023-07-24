#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -n 64 # Number of cores
#SBATCH --mem=256G
#SBATCH --job-name=run-example-experiment
#SBATCH -o run_example_experiment_slurm-%A_%a.out

python run_example_experiment.py
