"""
The purpose of this python3 script is to run a benchmark experiment on ACE runtime and memory.
"""


import datetime
import pandas as pd
import random
import tracemalloc
from acelib.main import run_ace_generate
from acelib.utilities import convert_dataframe_to_peptides


DESIGN_CSV_FILE = '/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/raw/runtime_memory_designs.csv'
OUTPUT_DIR = '/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/02_benchmark_ace_runtime_memory'
ITERATIONS = 10
GOLFY_MAX_ITERS = 2000
GOLFY_ALLOW_EXTRA_POOLS = False


def generate_peptides_dataframe(num_peptides: int):
    data = {
        'peptide_id': ['peptide_%i' % i for i in range(1, num_peptides + 1)],
        'peptide_sequence': ['' for _ in range(0, num_peptides)]
    }
    return pd.DataFrame(data)


if __name__ == '__main__':
    df_designs = pd.read_csv(DESIGN_CSV_FILE)
    data = {
        'num_peptides': [],
        'num_peptides_per_pool': [],
        'num_coverage': [],
        'iteration': [],
        'runtime_seconds': [],
        'peak_memory_mb': []
    }
    for index, value in df_designs.iterrows():
        num_peptides = value['num_peptides']
        num_peptides_per_pool = value['num_peptides_per_pool']
        num_coverage = value['num_coverage']
        df_peptides = generate_peptides_dataframe(num_peptides=num_peptides)
        peptides = convert_dataframe_to_peptides(df_peptides=df_peptides)

        for i in range(0, ITERATIONS):
            # Runtime
            start_datetime = datetime.datetime.now()
            run_ace_generate(
                peptides=peptides,
                num_peptides_per_pool=num_peptides_per_pool,
                num_coverage=num_coverage,
                cluster_peptides=False,
                trained_model_file='',
                mode='golfy',
                golfy_max_iters=GOLFY_MAX_ITERS,
                golfy_allow_extra_pools=GOLFY_ALLOW_EXTRA_POOLS,
                verbose=False
            )
            end_datetime = datetime.datetime.now()
            runtime_seconds = (end_datetime - start_datetime).total_seconds()

            # Memory
            tracemalloc.start()
            run_ace_generate(
                peptides=peptides,
                num_peptides_per_pool=num_peptides_per_pool,
                num_coverage=num_coverage,
                cluster_peptides=False,
                trained_model_file='',
                mode='golfy',
                golfy_max_iters=GOLFY_MAX_ITERS,
                golfy_allow_extra_pools=GOLFY_ALLOW_EXTRA_POOLS,
                verbose=False
            )
            peak_memory_mb = tracemalloc.get_traced_memory()[1] / (10 ** 6)
            tracemalloc.stop()

            data['num_peptides'].append(num_peptides)
            data['num_peptides_per_pool'].append(num_peptides_per_pool)
            data['num_coverage'].append(num_coverage)
            data['iteration'].append(i)
            data['runtime_seconds'].append(runtime_seconds)
            data['peak_memory_mb'].append(peak_memory_mb)

    random_id = random.randint(100000000, 999999999)
    df_benchmark = pd.DataFrame(data)
    df_benchmark.to_csv('%s/ace_runtime_memory_benchmark_experiment_results_%i.csv' % (OUTPUT_DIR, random_id),
                        index=False)
