# ACE Configurator for ELISpot Analysis

## 01. Installation

We developed a package to simulate *in silico* experiments for `ACE`:

```shell
pip install . --verbose
```

## 02. Repository Structure
```
├── README.md
├── examples
│   └── run_example_experiment.py               # acesim example
├── experiments                                 # experiment codes 
│   ├── 01_benchmark_ace_vs_naive_approach
│   ├── 02_benchmark_ace_features
│   ├── 03_benchmark_speed_and_memory
│   └── 04_generate_resource_configurations
├── pyproject.toml
├── src                                         # acesim source codes
│   └── acesim
├── test                                        # acesim unit tests
│   ├── data
│   │   └── iedb_mmer_all.csv
│   ├── data.py
│   └── test_experiment.py
└── unittest.sh
```

## 03. Experiments

* `01_benchmark_ace_vs_naive_approach` \
Benchmark `ACE` against a naive approach (i.e. random assignment). \
Evaluation metrics: number of total pools, sensitivity, specificity, precision.
  * With and without random experimental effects
    * `mu_immunogenic`
    * `mu_nonimmunogenic`
    * `dispersion_factor`
    * `method`
    * `alpha`
  * With and without sequence similarity feature (Alanine scanning)
<br/><br/>

* `02_benchmark_ace_features` \
Benchmark various `ACE` features. \
Evaluation metrics: number of total pools, sensitivity, specificity, precision.
  * `--mode`
  * `--golfy-max-iters`
  * `--num-peptides-per-pool`
  * `--num-coverage`
  * `--sequence-similarity-threshold`
<br/><br/>

* `03_benchmark_speed_and_memory` \
Benchmark `ACE` runtime and peak memory usage.
<br/><br/>

* `04_generate_resource_configurations` \
Generate all precomputed `ACE` ELISpot configurations.
<br/></br>
