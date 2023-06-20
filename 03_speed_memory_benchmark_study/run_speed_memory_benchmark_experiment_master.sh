#!/bin/bash

for i in {1..20}
do
  sbatch run_speed_memory_benchmark_experiment.sh
done

