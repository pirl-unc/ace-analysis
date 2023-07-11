#!/usr/bin/python3

"""
The purpose of this python3 script is to run one iteration of speed and memory benchmark study.

Last updated date: June 13, 2023
"""


import datetime
import random
import pandas as pd
import tracemalloc
from ortools.sat.python import cp_model
from acelib.main import run_ace_generate


CONFIGURATION_TSV_FILE = '/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/raw/speed_memory_benchmark_configuration.tsv'
OUTPUT_DIR = '/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/01_speed_memory_benchmark_study'
NUM_PROCESSES = 48


def generate_peptides_dataframe(num_peptides: int):
    data = {
        'peptide_id': ['peptide_%i' % i for i in range(1, num_peptides + 1)],
        'peptide_sequence': ['' for _ in range(0, num_peptides)]
    }
    return pd.DataFrame(data)


if __name__ == '__main__':
    df_study_config = pd.read_csv(CONFIGURATION_TSV_FILE, sep='\t')
    data = {
        'num_peptides': [],
        'num_peptides_per_pool': [],
        'num_coverage': [],
        'duration_in_seconds': [],
        'peak_memory_in_mb': [],
        'successful': []
    }
    benchmark_index = random.randint(1000000, 9999999)
    for index, value in df_study_config.iterrows():
        num_peptides = value['num_peptides']
        num_peptides_per_pool = value['num_peptides_per_pool']
        num_coverage = value['num_coverage']
        df_peptides = generate_peptides_dataframe(num_peptides=num_peptides)

        # Speed
        start_datetime = datetime.datetime.now()
        status, df_configuration = run_ace_generate(
            df_peptides=df_peptides,
            num_peptides_per_pool=num_peptides_per_pool,
            num_coverage=num_coverage,
            num_processes=NUM_PROCESSES
        )
        end_datetime = datetime.datetime.now()
        duration_in_seconds = (end_datetime - start_datetime).total_seconds()

        # Memory
        tracemalloc.start()
        run_ace_generate(
            df_peptides=df_peptides,
            num_peptides_per_pool=num_peptides_per_pool,
            num_coverage=num_coverage,
            num_processes=NUM_PROCESSES
        )
        mem = tracemalloc.get_traced_memory()[1] / (10 ** 6)
        tracemalloc.stop()

        data['num_peptides'].append(num_peptides)
        data['num_peptides_per_pool'].append(num_peptides_per_pool)
        data['num_coverage'].append(num_coverage)
        data['duration_in_seconds'].append(duration_in_seconds)
        data['peak_memory_in_mb'].append(mem)

        if status == cp_model.OPTIMAL:
            data['successful'].append(True)
            output_tsv_file = '%s/benchmark%i_%ipeptides_%iperpool_%ix.csv' % \
                              (OUTPUT_DIR, benchmark_index, num_peptides, num_peptides_per_pool, num_coverage)
            df_configuration.to_csv(output_tsv_file, index=False)
        else:
            data['successful'].append(False)

    df_benchmark = pd.DataFrame(data)
    output_tsv_file = '%s/benchmark%i_results.csv' % (OUTPUT_DIR, benchmark_index)
    df_benchmark.to_csv(output_tsv_file, index=False)
