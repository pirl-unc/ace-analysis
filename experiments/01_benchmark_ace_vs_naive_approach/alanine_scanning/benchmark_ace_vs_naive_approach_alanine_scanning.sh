#!/bin/bash

#SBATCH --time=72:00:00
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -n 64 # Number of cores
#SBATCH --mem=256G
#SBATCH --job-name=benchmark-ace-vs-naive-approach-alanine-scanning
#SBATCH -o benchmark_ace_vs_naive_approach_alanine_scanning_slurm-%A_%a.out

export OMP_NUM_THREADS=1
python benchmark_ace_vs_naive_approach_alanine_scanning.py
