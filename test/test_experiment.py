import random
import pandas as pd
from .data import get_data_path
from acesim.experiment import Experiment


def test_sample_peptides(sampled_peptides):
    print(sampled_peptides.head(n=5))


def test_simulate_peptide_spot_counts(sampled_peptides_spot_counts):
    print(sampled_peptides_spot_counts.head(n=5))


def test_aggregate_pool_spot_counts(pool_spot_counts):
   print(pool_spot_counts.head(n=5))


def test_identify_hit_pools_threshold(pool_spot_counts):
    return Experiment.identify_hit_pools(
        df_readouts=pool_spot_counts,
        num_peptides_per_pool=10,
        method='threshold',
        threshold=50.0,
        mu_nonimmunogenic=10.0,
        dispersion_factor=1.0,
        alpha=0.05,
        negative_control_sampling_iters=100
    )


def test_identify_hit_pools_adaptive(pool_spot_counts):
    return Experiment.identify_hit_pools(
        df_readouts=pool_spot_counts,
        num_peptides_per_pool=10,
        method='adaptive',
        threshold=50.0,
        mu_nonimmunogenic=10.0,
        dispersion_factor=1.0,
        alpha=0.05,
        negative_control_sampling_iters=100
    )


def test_identify_hit_pools_combined(pool_spot_counts):
    return Experiment.identify_hit_pools(
        df_readouts=pool_spot_counts,
        num_peptides_per_pool=10,
        method='combined',
        threshold=50.0,
        mu_nonimmunogenic=10.0,
        dispersion_factor=1.0,
        alpha=0.05,
        negative_control_sampling_iters=100
    )

