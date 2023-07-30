#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -n 2 # Number of cores
#SBATCH --mem=64G
#SBATCH --job-name=benchmark-ace-vs-naive-approach-deconvolutethis-designs-results-aggregation
#SBATCH -o benchmark_ace_vs_naive_approach_deconvolutethis_designs_results_aggregation_slurm.out

python benchmark_ace_vs_naive_approach_deconvolutethis_designs_results_aggregation.py
