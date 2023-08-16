#!/bin/bash
for i in {1..10}
do
   sbatch 01_run_reduced_designs_experiments.sh
done