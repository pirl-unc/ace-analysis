import random
import pandas as pd
import pytest
from .data import get_data_path
from acesim.experiment import Experiment
from acesim.configurator_ace import AceConfigurator


@pytest.fixture
def sampled_peptides() -> pd.DataFrame:
    df_ref_peptides = pd.read_csv(get_data_path(name='held_out_data_w_negatives.csv'))
    df_peptides = Experiment.sample_peptides(
        df_ref_peptides=df_ref_peptides,
        num_peptides=100,
        num_peptides_immunogenic=10,
        peptide_length=9,
        sampling_method='random'
    )
    return df_peptides


@pytest.fixture
def sampled_peptides_spot_counts(sampled_peptides) -> pd.DataFrame:
    df_peptides_spot_counts = Experiment.simulate_peptide_spot_counts(
        df_peptides=sampled_peptides,
        num_coverage=3,
        random_effects=True,
        mu_immunogenic=100.0,
        mu_nonimmunogenic=10.0,
        dispersion_factor=1.0
    )
    return df_peptides_spot_counts


@pytest.fixture
def pool_spot_counts(sampled_peptides_spot_counts):
    amino_acids = [
        "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
        "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"
    ]
    peptides = []
    for i in range(0, 100):
        sampled_amino_acids = random.sample(amino_acids, 9)
        peptides.append(('peptide_%i' % (i+1), ''.join(sampled_amino_acids)))
    ace_configurator = AceConfigurator(name='')
    block_assignment, preferred_neighbors = ace_configurator.generate_assignment(
        peptides=peptides,
        num_peptides_per_pool=10,
        num_coverage=3,
        trained_model_file='',
        cluster_peptides=False,
        mode='golfy',
        sequence_similarity_threshold=0.0,
        sequence_similarity_function='',
        golfy_random_seed=1,
        golfy_strategy='random',
        golfy_max_iters=1000,
        golfy_allow_extra_pools=False,
        plate_size=96
    )
    df_readout = Experiment.aggregate_pool_spot_counts(
        df_assignments=block_assignment.to_dataframe(),
        df_peptides_spot_counts=sampled_peptides_spot_counts,
        mu_nonimmunogenic=10.0,
        dispersion_factor=1.0,
        false_negative_rate=0.0,
        num_peptides_per_pool=10
    )
    return df_readout
