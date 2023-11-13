#!/bin/bash

#SBATCH --time=72:00:00
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -n 64 # Number of cores
#SBATCH --mem=256G
#SBATCH --job-name=run-tiling-window-covid-spike-9mers-experiments
#SBATCH -o 03_run_tiling_window_covid_spike_9mers_experiments_slurm-%A_%a.out

export OMP_NUM_THREADS=1
python 03_run_tiling_window_covid_spike_9mers_experiments.py
