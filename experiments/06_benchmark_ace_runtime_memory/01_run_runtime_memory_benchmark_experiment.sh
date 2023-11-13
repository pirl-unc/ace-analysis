#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -n 64 # Number of cores
#SBATCH --mem=256G
#SBATCH --job-name=run-ace-runtime-memory-benchmark-experiment
#SBATCH -o 01_run_runtime_memory_benchmark_experiment_slurm.out

python 01_run_runtime_memory_benchmark_experiment.py
