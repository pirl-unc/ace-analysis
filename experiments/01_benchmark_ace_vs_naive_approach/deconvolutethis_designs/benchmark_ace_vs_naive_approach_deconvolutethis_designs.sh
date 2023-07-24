#!/bin/bash

#SBATCH --time=96:00:00
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -n 64 # Number of cores
#SBATCH --mem=256G
#SBATCH --job-name=benchmark-ace-vs-naive-approach-deconvolutethis-designs
#SBATCH -o benchmark_ace_vs_naive_approach_deconvolutethis_designs_slurm-%A_%a.out

export OMP_NUM_THREADS=1
python benchmark_ace_vs_naive_approach_deconvolutethis_designs.py
