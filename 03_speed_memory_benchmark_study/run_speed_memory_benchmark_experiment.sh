#!/bin/bash

#SBATCH --time=120:00:00
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -n 50 # Number of cores
#SBATCH --mem=256G
#SBATCH --job-name=run-ace-speed-memory-benchmark-experiment
#SBATCH -o run_speed_memory_benchmark_experiment_slurm.out

python run_speed_memory_benchmark_experiment.py
