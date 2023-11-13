#!/bin/bash

#SBATCH --time=2:00:00
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -n 2 # Number of cores
#SBATCH --mem=16G
#SBATCH --job-name=aggregate-tiling-window-covid-spike-9mers-experiments-results
#SBATCH -o 04_aggregate_tiling_window_covid_spike_9mers_experiments_results_slurm.out

python 04_aggregate_tiling_window_covid_spike_9mers_experiments_results.py
