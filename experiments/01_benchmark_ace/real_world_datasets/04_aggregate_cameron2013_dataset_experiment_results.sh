#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -n 2 # Number of cores
#SBATCH --mem=64G
#SBATCH --job-name=aggregate-cameron2013-dataset-experiment-results
#SBATCH -o 04_aggregate_cameron2013_dataset_experiment_results_slurm.out

python 04_aggregate_cameron2013_dataset_experiment_results.py
