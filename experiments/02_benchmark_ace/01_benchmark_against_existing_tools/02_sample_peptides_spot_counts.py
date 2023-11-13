"""
The purpose of this python3 script is to sample spot counts for peptides.
"""


import glob
import os
import pandas as pd
import multiprocessing as mp
from acesim.experiment import Experiment


NUM_PROCESSES = 64
PEPTIDE_SAMPLING_COVERAGE_SCHEME_TSV_FILE = "/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/raw/benchmark_against_existing_tools/peptide_sampling_coverage_scheme.tsv"
SAMPLED_PEPTIDES_DIR = "/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/02_benchmark_ace/01_benchmark_against_existing_tools/01_sampled_peptides"
OUTPUT_DIR = "/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/02_benchmark_ace/01_benchmark_against_existing_tools/02_sampled_peptides_spot_counts"


def worker(task):
    df_scheme = task[0]
    tsv_file = task[1]
    mu_immunogenic = task[2]
    mu_nonimmunogenic = task[3]
    dispersion_factor = task[4]
    output_dir = task[5]
    tsv_file_basename = os.path.basename(tsv_file)
    df_peptides = pd.read_csv(tsv_file, sep='\t')
    num_peptides = len(df_peptides)
    df_matched = df_scheme.loc[df_scheme['num_peptides'] == num_peptides, :]
    for _, row in df_matched.iterrows():
        num_coverage = row['num_coverage']
        df_peptides_spot_counts = Experiment.simulate_peptide_spot_counts(
            df_peptides=df_peptides,
            num_coverage=num_coverage,
            random_effects=True,
            mu_immunogenic=mu_immunogenic,
            mu_nonimmunogenic=mu_nonimmunogenic,
            dispersion_factor=dispersion_factor
        )
        output_file = output_dir + '/' + \
                      tsv_file_basename.replace('.tsv','') + '_' + \
                      str(num_coverage) + 'x_spot_counts.tsv'
        df_peptides_spot_counts.to_csv(output_file, sep='\t', index=False)


if __name__ == "__main__":
    print("Started")
    if not os.path.exists(OUTPUT_DIR + '/300immunogenic_30nonimmunogenic_1dispersion/'):
        os.makedirs(OUTPUT_DIR + '/300immunogenic_30nonimmunogenic_1dispersion/')
    if not os.path.exists(OUTPUT_DIR + '/300immunogenic_30nonimmunogenic_0dispersion/'):
        os.makedirs(OUTPUT_DIR + '/300immunogenic_30nonimmunogenic_0dispersion/')
    pool = mp.Pool(processes=NUM_PROCESSES)
    tasks = []
    df_scheme = pd.read_csv(PEPTIDE_SAMPLING_COVERAGE_SCHEME_TSV_FILE, sep='\t')
    for tsv_file in glob.glob(SAMPLED_PEPTIDES_DIR + '/*.tsv'):
        tasks.append((df_scheme, tsv_file,  300.0, 30.0, 1.0, OUTPUT_DIR + '/300immunogenic_30nonimmunogenic_1dispersion'))
        tasks.append((df_scheme, tsv_file,  300.0, 30.0, 0.0, OUTPUT_DIR + '/300immunogenic_30nonimmunogenic_0dispersion'))
    pool.map(worker, tasks)
    pool.close()
    pool.join()
    print("Finished successfully")
