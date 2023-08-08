#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -n 10 # Number of cores
#SBATCH --mem=40G
#SBATCH --job-name=run-design-exploration-experiment
#SBATCH -o 01_run_design_exploration_experiment_slurm-%A_%a.out

export OMP_NUM_THREADS=1
python 01_run_design_exploration_experiment.py
