"""
The purpose of this python3 script is to sample peptides.
"""


import os
import pandas as pd
import multiprocessing as mp
from acesim.experiment import Experiment


NUM_PROCESSES = 32
PEPTIDE_SAMPLING_SCHEME_TSV_FILE = "/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/raw/benchmark_higher_positivity_rates/peptide_sampling_scheme.tsv"
IEDB_CSV_FILE = "/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/raw/references/large_dummy_sequence_dataset.csv"
OUTPUT_DIR = "/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/02_benchmark_ace/02_benchmark_higher_positivity_rates/01_sampled_peptides"
PEPTIDE_LENGTH = 9
SAMPLING_METHOD = 'levenshtein'
NUM_REPLICATES = 100


def worker(task):
    df_ref_peptides = task[0]
    num_peptides = task[1]
    num_peptides_immunogenic = task[2]
    replicate = task[3]
    df_peptides = Experiment.sample_peptides(
        df_ref_peptides=df_ref_peptides,
        num_peptides=num_peptides,
        num_peptides_immunogenic=num_peptides_immunogenic,
        peptide_length=PEPTIDE_LENGTH,
        sampling_method=SAMPLING_METHOD
    )
    df_peptides.to_csv('%s/%ipeptides_%iimmunogenic_rep%i.tsv' %
                       (OUTPUT_DIR, num_peptides, num_peptides_immunogenic, replicate),
                       sep='\t', index=False)


if __name__ == "__main__":
    print("Started")
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    df_ref_peptides = pd.read_csv(IEDB_CSV_FILE)
    df_scheme = pd.read_csv(PEPTIDE_SAMPLING_SCHEME_TSV_FILE, sep='\t')
    pool = mp.Pool(processes=NUM_PROCESSES)
    tasks = []
    for _, row in df_scheme.iterrows():
        num_peptides = int(row['num_peptides'])
        num_peptides_immunogenic = int(row['num_peptides_immunogenic'])
        for i in range(0, NUM_REPLICATES):
            tasks.append((df_ref_peptides, num_peptides, num_peptides_immunogenic, i+1))
    pool.map(worker, tasks)
    pool.close()
    pool.join()
    print("Finished successfully")
