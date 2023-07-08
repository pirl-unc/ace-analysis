#!/usr/bin/python3

"""
The purpose of this python3 script is to generate all ACE resource configurations.

Last updated date: July 5, 2023
"""


import pandas as pd
from acelib.main import run_ace_generate


CONFIGURATION_TSV_FILE = '/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/raw/resource_configurations.tsv'
OUTPUT_DIR = '/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/01_resource_configurations'
NUM_PROCESSES = 48
GOLFY_MAX_ITERS = 100000
RANDOM_SEED = 42
NUM_PEPTIDES_PER_BATCH = 100


def generate_peptides_dataframe(num_peptides: int):
    data = {
        'peptide_id': ['peptide_%i' % i for i in range(1, num_peptides + 1)],
        'peptide_sequence': ['' for _ in range(0, num_peptides)]
    }
    return pd.DataFrame(data)


if __name__ == '__main__':
    df_configurations = pd.read_csv(CONFIGURATION_TSV_FILE, sep='\t')
    data = {
        'num_peptides'
    }
    optimality = []
    for index, value in df_configurations.iterrows():
        num_peptides = value['num_peptides']
        num_peptides_per_pool = value['num_peptides_per_pool']
        num_coverage = value['num_coverage']
        df_peptides = generate_peptides_dataframe(num_peptides=num_peptides)
        status, df_configuration = run_ace_generate(
            df_peptides=df_peptides,
            num_peptides_per_pool=num_peptides_per_pool,
            num_coverage=num_coverage,
            num_processes=NUM_PROCESSES,
            random_seed=RANDOM_SEED,
            num_peptides_per_batch=NUM_PEPTIDES_PER_BATCH,
            golfy_max_iters=GOLFY_MAX_ITERS
        )
        output_csv_file = '%s/%ipeptides_%iperpool_%ix.csv' % (OUTPUT_DIR, num_peptides, num_peptides_per_pool, num_coverage)
        df_configuration.to_csv(output_csv_file, index=False)
        optimality.append(status)

    df_configurations['optimality'] = optimality
    df_configurations.to_csv(OUTPUT_DIR + '/configurations_optimality.csv', index=False)

