#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -n 64 # Number of cores
#SBATCH --mem=256G
#SBATCH --job-name=run-deconvolution-methods-experiments
#SBATCH -o 01_run_deconvolution_methods_experiments.out

export OMP_NUM_THREADS=1
python 01_run_deconvolution_methods_experiments.py
