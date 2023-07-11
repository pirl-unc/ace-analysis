#!/usr/bin/python3

"""
The purpose of this python3 script is to generate all ACE resource configurations.

Last updated date: July 11, 2023
"""


import pandas as pd
from acelib.main import run_ace_golfy, run_ace_sat_solver


CONFIGURATION_CSV_FILE = '/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/raw/resource_configurations.csv'
OUTPUT_DIR = '/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/04_generate_resource_configurations'
NUM_PROCESSES = 64
GOLFY_MAX_ITERS = 2000
GOLFY_INIT_MODE = 'greedy'
RANDOM_SEED = 42
NUM_PEPTIDES_PER_BATCH = 100


def generate_peptides_dataframe(num_peptides: int):
    data = {
        'peptide_id': ['peptide_%i' % i for i in range(1, num_peptides + 1)],
        'peptide_sequence': ['' for _ in range(0, num_peptides)]
    }
    return pd.DataFrame(data)


if __name__ == '__main__':
    df_configurations = pd.read_csv(CONFIGURATION_CSV_FILE)
    data = {
        'num_peptides': [],
        'num_peptides_per_pool': [],
        'num_coverage': [],
        'optimal': []
    }
    for index, value in df_configurations.iterrows():
        num_peptides = value['num_peptides']
        num_peptides_per_pool = value['num_peptides_per_pool']
        num_coverage = value['num_coverage']
        df_peptides = generate_peptides_dataframe(num_peptides=num_peptides)
        if num_coverage > 3 or num_peptides_per_pool > 10:
            is_valid, df_configuration = run_ace_golfy(
                df_peptides=df_peptides,
                num_peptides_per_pool=num_peptides_per_pool,
                num_coverage=num_coverage,
                random_seed=RANDOM_SEED,
                max_iters=GOLFY_MAX_ITERS,
                init_mode=GOLFY_INIT_MODE
            )
            data['num_peptides'].append(num_peptides)
            data['num_peptides_per_pool'].append(num_peptides_per_pool)
            data['num_coverage'].append(num_coverage)
            data['optimal'].append(is_valid)
        else:
            df_configuration = run_ace_sat_solver(
                df_peptides=df_peptides,
                num_peptides_per_pool=num_peptides_per_pool,
                num_coverage=num_coverage,
                num_peptides_per_batch=NUM_PEPTIDES_PER_BATCH,
                random_seed=RANDOM_SEED,
                num_processes=NUM_PROCESSES
            )
            data['num_peptides'].append(num_peptides)
            data['num_peptides_per_pool'].append(num_peptides_per_pool)
            data['num_coverage'].append(num_coverage)
            data['optimal'].append(True)

        output_csv_file = '%s/%ipeptides_%iperpool_%ix.csv' % (OUTPUT_DIR, num_peptides, num_peptides_per_pool, num_coverage)
        df_configuration.to_csv(output_csv_file, index=False)

    df_run = pd.DataFrame(data)
    df_run.to_csv(OUTPUT_DIR + '/configurations_optimality.csv', index=False)

