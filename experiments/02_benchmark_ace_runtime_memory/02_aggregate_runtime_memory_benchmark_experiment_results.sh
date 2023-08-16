#!/bin/bash

#SBATCH --time=1:00:00
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -n 2 # Number of cores
#SBATCH --mem=2G
#SBATCH --job-name=aggregate-runtime-memory-experiment-results
#SBATCH -o 02_aggregate_runtime_memory_benchmark_experiment_results_slurm.out

python 02_aggregate_runtime_memory_benchmark_experiment_results.py
