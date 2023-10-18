#!/bin/bash

#SBATCH --time=72:00:00
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -n 1 # Number of cores
#SBATCH --mem=32G
#SBATCH --job-name=run-ace-runtime-memory-benchmark-experiment
#SBATCH -o 01_run_runtime_memory_benchmark_experiment_slurm.out

python 01_run_runtime_memory_benchmark_experiment.py
