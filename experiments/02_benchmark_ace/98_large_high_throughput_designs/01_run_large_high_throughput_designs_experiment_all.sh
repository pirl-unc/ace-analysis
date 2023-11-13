#!/bin/bash
for i in {1..10}
do
   sbatch 01_run_large_high_throughput_designs_experiment.sh
done