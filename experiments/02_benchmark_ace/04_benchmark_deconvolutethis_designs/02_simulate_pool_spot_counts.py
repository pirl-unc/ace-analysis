"""
The purpose of this python3 script is to simulate pool spot counts.
"""


import glob
import openpyxl
import os
import pandas as pd
import numpy as np
import multiprocessing as mp
from acesim.experiment import Experiment


NUM_PROCESSES = 64
DESIGN_CONFIGURATIONS_TSV_FILE = "/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/raw/benchmark_deconvolutethis_designs/design_configurations.tsv"
SAMPLED_PEPTIDES_SPOT_COUNTS_DIR = "/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/02_benchmark_ace/01_benchmark_against_existing_tools/02_sampled_peptides_spot_counts/300immunogenic_30nonimmunogenic_0dispersion"   # reuse sampled peptides
CONFIGURATIONS_DIR = "/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/02_benchmark_ace/04_benchmark_deconvolutethis_designs/01_generated_configurations"
OUTPUT_DIR = "/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/02_benchmark_ace/04_benchmark_deconvolutethis_designs/02_simulated_pool_spot_counts/300immunogenic_30nonimmunogenic_0dispersion_fnr0.00"
MU_IMMUNOGENIC = 300.0
MU_NONIMMUNOGENIC = 30.0
DISPERSION_FACTOR = 0.0
SATURATED_SPOT_COUNT = 600
FALSE_NEGATIVE_RATE = 0.00
NUM_NEGATIVE_CONTROL_WELLS = 3


def simulate_pool_spot_counts_worker(task):
    design_tsv_file = task[0]
    peptides_spot_counts_tsv_file = task[1]
    output_file = task[2]
    num_peptides_per_pool = task[3]
    df_assignments = pd.read_csv(design_tsv_file, sep='\t')
    df_peptides_spot_counts = pd.read_csv(peptides_spot_counts_tsv_file, sep='\t')
    df_readout = Experiment.aggregate_pool_spot_counts(
        df_assignments=df_assignments,
        df_peptides_spot_counts=df_peptides_spot_counts,
        mu_nonimmunogenic=MU_NONIMMUNOGENIC,
        dispersion_factor=DISPERSION_FACTOR,
        false_negative_rate=FALSE_NEGATIVE_RATE,
        num_peptides_per_pool=num_peptides_per_pool
    )
    data = {
        'pool_id': [],
        'plate_id': [],
        'well_id': [],
        'spot_count': []
    }
    for _, row in df_readout.iterrows():
        df_matched = df_assignments.loc[df_assignments['pool_id'] == row['pool_id'],:]
        data['pool_id'].append(int(row['pool_id']))
        data['spot_count'].append(round(row['spot_count']))
        data['plate_id'].append(df_matched['plate_id'].values.tolist()[0])
        data['well_id'].append(df_matched['well_id'].values.tolist()[0])

    # Simulate negative control well spot counts
    for _ in range(0, NUM_NEGATIVE_CONTROL_WELLS):
        negative_pool_spot_count = Experiment.sample_spot_counts(
            mean=MU_NONIMMUNOGENIC,
            dispersion_factor=DISPERSION_FACTOR,
            num_coverage=num_peptides_per_pool
        )
        data['pool_id'].append('negative')
        data['plate_id'].append(1)
        data['well_id'].append('')
        data['spot_count'].append(round(np.sum(np.array(negative_pool_spot_count) / float(num_peptides_per_pool))))
    df_readout_ = pd.DataFrame(data)
    df_readout_.to_csv(output_file, sep='\t', index=False)


def simulate_pool_spot_counts():
    # Step 1. Create tasks
    df_designs = pd.read_csv(DESIGN_CONFIGURATIONS_TSV_FILE, sep='\t')
    tasks = []
    for _, row in df_designs.iterrows():
        num_peptides = row['num_peptides']
        num_peptides_per_pool = row['num_peptides_per_pool']
        num_coverage = row['num_coverage']
        num_immunogenic_peptides = row['num_immunogenic_peptides']

        output_dir_ = '%s/ace-s/%ipeptides_%iperpool_%ix' % (OUTPUT_DIR, num_peptides, num_peptides_per_pool, num_coverage)
        if not os.path.exists(output_dir_):
            os.makedirs(output_dir_)

        # ACE-S design
        design_tsv_files = glob.glob('%s/ace-s/%ipeptides_%iperpool_%ix/%ipeptides_%iimmunogenic_rep*_%iperpool_%ix_ace-s_design_*.tsv' %
                                     (CONFIGURATIONS_DIR,
                                      num_peptides, num_peptides_per_pool, num_coverage,
                                      num_peptides, num_immunogenic_peptides, num_peptides_per_pool, num_coverage))
        for design_tsv_file in design_tsv_files:
            design_tsv_file_basename = os.path.basename(design_tsv_file)
            design_tsv_file_basename_elements = design_tsv_file_basename.split('_')
            replicate = design_tsv_file_basename_elements[2]
            peptides_spots_count_file = '%s/%ipeptides_%iimmunogenic_%s_%ix_spot_counts.tsv' % (SAMPLED_PEPTIDES_SPOT_COUNTS_DIR, num_peptides, num_immunogenic_peptides, replicate, num_coverage)
            output_file = ('%s/%s' % (output_dir_, design_tsv_file_basename.replace('.tsv', '_pool_spot_counts.tsv')))
            tasks.append((design_tsv_file, peptides_spots_count_file, output_file, num_peptides_per_pool))

    # Step 2. Multiprocess tasks
    pool = mp.Pool(processes=NUM_PROCESSES)
    pool.map(simulate_pool_spot_counts_worker, tasks)
    pool.close()
    pool.join()


if __name__ == "__main__":
    print("Started")
    simulate_pool_spot_counts()
    print("Finished successfully")
