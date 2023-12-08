# ACE Configurator for ELISpot Analysis

## 01. Installation

We developed a package to simulate *in silico* experiments for `ACE`:

```shell
pip install . --verbose
```

## 02. Repository Structure
```
├── README.md
├── experiments                                   # experiment codes 
│   ├── 01_precision_recall_violations_iterations
│   ├── 02_benchmark_ace
│   ├── 03_alanine_scanning
│   ├── 04_real_world_datasets
│   ├── 05_design_exploration
│   └── 06_benchmark_ace_runtime_memory
├── figures                                       # codes used to generate figures
├── src                                           # acesim source codes
│   └── acesim
├── tables                                        # codes used to generate tables
├── test                                          # acesim unit tests
└── unittest.sh
```
