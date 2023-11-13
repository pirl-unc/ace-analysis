"""
The purpose of this python3 script is to run a benchmark experiment on ACE runtime and memory.
"""


import datetime
import pandas as pd
import random
import tracemalloc
import multiprocessing as mp
import os
from acelib.main import run_ace_generate


REFERENCE_PEPTIDES_CSV_FILE = '/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/raw/runtime_memory/strict_holdout_negatives.csv'
DESIGN_CSV_FILE = '/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/raw/runtime_memory/runtime_memory_designs.csv'
TRAINED_MODEL_FILE = "/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/raw/trained_models/trained_model_w_data_augmentation_b3000.pt"
OUTPUT_DIR = '/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/06_benchmark_ace_runtime_memory'
NUM_PROCESSES = 64
NUM_ITERATIONS = 100
GOLFY_MAX_ITERS = 2000
GOLFY_ALLOW_EXTRA_POOLS = False


def worker(task):
    peptides = task[0]
    num_peptides = task[1]
    num_peptides_per_pool = task[2]
    num_coverage = task[3]
    rep_id = task[4]
    data = {
        'num_peptides': [],
        'num_peptides_per_pool': [],
        'num_coverage': [],
        'iteration': [],
        'runtime_seconds': [],
        'peak_memory_mb': []
    }
    start_datetime = datetime.datetime.now()
    tracemalloc.start()
    run_ace_generate(
        peptides=peptides,
        num_peptides_per_pool=num_peptides_per_pool,
        num_coverage=num_coverage,
        cluster_peptides=True,
        trained_model_file=TRAINED_MODEL_FILE,
        mode='golfy',
        golfy_max_iters=GOLFY_MAX_ITERS,
        golfy_allow_extra_pools=GOLFY_ALLOW_EXTRA_POOLS,
        verbose=False
    )
    peak_memory_mb = tracemalloc.get_traced_memory()[1] / (10 ** 6)
    tracemalloc.stop()
    end_datetime = datetime.datetime.now()
    runtime_seconds = (end_datetime - start_datetime).total_seconds()

    data['num_peptides'].append(num_peptides)
    data['num_peptides_per_pool'].append(num_peptides_per_pool)
    data['num_coverage'].append(num_coverage)
    data['iteration'].append(rep_id)
    data['runtime_seconds'].append(runtime_seconds)
    data['peak_memory_mb'].append(peak_memory_mb)

    return pd.DataFrame(data)


if __name__ == '__main__':
    print("Started")
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    df_ref_peptides = pd.read_csv(REFERENCE_PEPTIDES_CSV_FILE)
    df_designs = pd.read_csv(DESIGN_CSV_FILE)
    tasks = []
    for index, value in df_designs.iterrows():
        num_peptides = value['num_peptides']
        num_peptides_per_pool = value['num_peptides_per_pool']
        num_coverage = value['num_coverage']
        for rep_id in range(0, NUM_ITERATIONS):
            df_peptides = df_ref_peptides.sample(n=num_peptides)
            peptide_idx = 1
            peptides = []
            for _, row in df_peptides.iterrows():
                peptides.append(('peptide_%i' % peptide_idx, row['epitope']))
                peptide_idx += 1
            tasks.append((peptides, num_peptides, num_peptides_per_pool, num_coverage, rep_id + 1))
    pool = mp.Pool(processes=NUM_PROCESSES)
    results = pool.map(worker, tasks)
    pool.close()
    pool.join()

    df_results = pd.DataFrame()
    for result in results:
        df_results = pd.concat([df_results, result], axis=0, ignore_index=True)
    df_results.to_csv(OUTPUT_DIR + '/runtime_memory_experiments_results.tsv', sep='\t', index=False)
    print("Finished successfully")
