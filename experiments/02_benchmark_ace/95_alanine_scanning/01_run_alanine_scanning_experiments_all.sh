#!/bin/bash
for i in {1..10}
do
   sbatch 01_run_alanine_scanning_experiments.sh
done