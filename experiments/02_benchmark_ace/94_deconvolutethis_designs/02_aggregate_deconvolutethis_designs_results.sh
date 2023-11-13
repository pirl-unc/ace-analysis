#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -n 2 # Number of cores
#SBATCH --mem=64G
#SBATCH --job-name=aggregate-deconvolutethis-designs-experiment-results
#SBATCH -o 02_aggregate_deconvolutethis_designs_results_slurm.out

python 02_aggregate_deconvolutethis_designs_results.py
