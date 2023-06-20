#!/bin/bash

#SBATCH --time=72:00:00
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -n 50 # Number of cores
#SBATCH --mem=256G
#SBATCH --job-name=investigate-hit-rate-vs-candidate-hits
#SBATCH -o investigate_hit_rate_vs_candidate_hits_slurm.out

python investigate_hit_rate_vs_candidate_hits.py
