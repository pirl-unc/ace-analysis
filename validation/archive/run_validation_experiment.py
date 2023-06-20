#!/usr/bin/python3

"""
The purpose of this python3 script is to run a validation study of ACE configurations.

Author: Jin Seok (Andy) Lee, Dhuvarakesh Karthikeyan

Last updated date: Aug 15, 2022
"""


import pandas as pd
from acelib.logger import *
from acelib.solver import *
from acelib.verification import *
from acelib.logger import *
from acelib.visualization import *


VALIDATION_STUDY_CONFIGURATION_TSV_FILE = "/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/raw/validation_study_configurations.tsv"
OUTPUT_DIR = "/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/validation_experiments"
NUM_CORES = 64


if __name__ == '__main__':
    logger = get_logger(__name__)

    # Step 1. Read the validation study configuration
    df_study_config = pd.read_csv(VALIDATION_STUDY_CONFIGURATION_TSV_FILE, sep='\t')

    # Step 2. Iterate through each configuration and run ACE
    logger.info("Started running validation experiment")
    for index, row in df_study_config.iterrows():
        n_peptides = int(row['n_peptides'])
        n_positive_peptides = int(row['n_positive_peptides'])
        n_coverage = int(row['n_coverage'])
        n_peptides_per_pool = int(row['n_peptides_per_pool'])
        n_simulations = int(row['n_simulations'])
        logger.info("%i peptides. %i positive peptides. %i coverage. %i peptides per pool." % (n_peptides, n_positive_peptides, n_coverage, n_peptides_per_pool))

        # for i in range(0, n_simulations):
        df_configuration = generate_assay_configuration(
            n_peptides=n_peptides,
            n_peptides_per_pool=n_peptides_per_pool,
            n_coverage=n_coverage,
            num_threads=NUM_CORES
        )
        output_file_prefix = OUTPUT_DIR + '/' + \
                             str(n_peptides) + '_peptides_' + \
                             str(n_peptides_per_pool) + '_peptides_per_pool_' + \
                             str(n_coverage) + 'x'
        verify_configuration_constraints(df_configuration=df_configuration,
                                         n_peptides=n_peptides,
                                         n_peptides_per_pool=n_peptides_per_pool,
                                         n_coverage=n_coverage)
        df_configuration.to_csv(output_file_prefix + '.tsv', sep='\t', index=False)
        plot_configuration_table(df_configuration=df_configuration,
                                 save_figure=True,
                                 output_pdf_file=output_file_prefix + '.pdf')
        print("\n")

    logger.info("Finished running validation experiment")

