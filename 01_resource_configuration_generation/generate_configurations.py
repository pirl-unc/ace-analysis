#!/usr/bin/python3

"""
The purpose of this python3 script is to generate all ACE resource configurations.

Last updated date: June 19, 2023
"""


import pandas as pd
from ortools.sat.python import cp_model
from acelib.main import run_ace_generate


CONFIGURATION_TSV_FILE = '/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/raw/ace_configurations.tsv'
OUTPUT_DIR = '/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/01_resource_configuration_generation'
NUM_PROCESSES = 48


def generate_peptides_dataframe(num_peptides: int):
    data = {
        'peptide_id': ['peptide_%i' % i for i in range(1, num_peptides + 1)],
        'peptide_sequence': ['' for _ in range(0, num_peptides)]
    }
    return pd.DataFrame(data)


if __name__ == '__main__':
    df_configurations = pd.read_csv(CONFIGURATION_TSV_FILE, sep='\t')
    for index, value in df_configurations.iterrows():
        num_peptides = value['num_peptides']
        num_peptides_per_pool = value['num_peptides_per_pool']
        num_coverage = value['num_coverage']
        df_peptides = generate_peptides_dataframe(num_peptides=num_peptides)

        # Speed
        status, df_configuration = run_ace_generate(
            df_peptides=df_peptides,
            num_peptides_per_pool=num_peptides_per_pool,
            num_coverage=num_coverage,
            num_processes=NUM_PROCESSES
        )

        if status == cp_model.OPTIMAL:
            output_tsv_file = '%s/%ipeptides_%iperpool_%ix.csv' % \
                              (OUTPUT_DIR, num_peptides, num_peptides_per_pool, num_coverage)
            df_configuration.to_csv(output_tsv_file, index=False)
        else:
            print("Optimal configuration could not be found for the following configuration: %i peptides, %i peptides per pool, %ix coverage." %
                  (num_peptides, num_peptides_per_pool, num_coverage))

