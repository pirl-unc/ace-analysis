"""
The purpose of this python3 script is to deconvolve pool spot counts.
"""


import glob
import os.path
import pandas as pd
import numpy as np
import multiprocessing as mp
from acelib.block_assignment import BlockAssignment
from acesim.deconvolver_ace import AceDeconvolver


NUM_PROCESSES = 64
CONFIGURATIONS_DIR = "/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/02_benchmark_ace/04_benchmark_deconvolutethis_designs/01_generated_configurations"
POOL_SPOT_COUNTS_DIR = "/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/02_benchmark_ace/04_benchmark_deconvolutethis_designs/02_simulated_pool_spot_counts/300immunogenic_30nonimmunogenic_0dispersion_fnr0.00"
OUTPUT_DIR = "/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/02_benchmark_ace/04_benchmark_deconvolutethis_designs/03_deconvolution_results/300immunogenic_30nonimmunogenic_0dispersion_fnr0.00"
NEGATIVE_CONTROL_MULTIPLIER = 3


def run_ace_deconvolution(
        pool_spot_counts_tsv_file: str,
        design_tsv_file: str,
        method: str,
        min_coverage: int,
        output_file: str
):
    df_readout = pd.read_csv(pool_spot_counts_tsv_file, sep='\t')
    min_pool_spot_count = np.mean(df_readout.loc[df_readout['pool_id'] == 'negative','spot_count'].values.tolist()) * NEGATIVE_CONTROL_MULTIPLIER
    df_readout = df_readout.loc[df_readout['pool_id'] != 'negative',:]
    df_readout['pool_id'] = df_readout['pool_id'].astype('int')
    df_assignments = pd.read_csv(design_tsv_file, sep='\t')
    df_assignments['pool_id'] = df_assignments['pool_id'].astype('int')
    block_assignment = BlockAssignment.load_from_dataframe(df_assignments=df_assignments)
    ace_deconvolver = AceDeconvolver(name='ace')
    deconvolution_result = ace_deconvolver.deconvolve(
        df_readout=df_readout,
        block_assignment=block_assignment,
        method=method,
        min_coverage=min_coverage,
        min_pool_spot_count=min_pool_spot_count
    )
    deconvolution_result.to_dataframe().to_csv(output_file, sep='\t', index=False)


if __name__ == "__main__":
    print("Started")

    pool = mp.Pool(processes=NUM_PROCESSES)

    # Step 1. ACE-S configurations - ACE-CEM
    for pool_spot_counts_tsv_file in glob.glob(POOL_SPOT_COUNTS_DIR + '/ace-s/*/*pool_spot_counts.tsv'):
        pool_spot_counts_tsv_file_basename = os.path.basename(pool_spot_counts_tsv_file)
        num_peptides = int(pool_spot_counts_tsv_file_basename.split('_')[0].replace('peptides',''))
        num_peptides_per_pool = int(pool_spot_counts_tsv_file_basename.split('_')[3].replace('perpool',''))
        num_coverage = int(pool_spot_counts_tsv_file_basename.split('_')[4].replace('x',''))
        config_dir_ = '%s/ace-s/%ipeptides_%iperpool_%ix' % (CONFIGURATIONS_DIR, num_peptides, num_peptides_per_pool, num_coverage)
        design_tsv_file = glob.glob(config_dir_ + '/' + '_'.join(pool_spot_counts_tsv_file_basename.split('_')[0:7]) + '*.tsv')[0]
        output_dir_ = '%s/%ipeptides_%iperpool_%ix' % (OUTPUT_DIR, num_peptides, num_peptides_per_pool, num_coverage)
        if not os.path.exists(output_dir_ + '/ace-s/ace_cem/'):
            os.makedirs(output_dir_ + '/ace-s/ace_cem/')
        output_file = output_dir_ + '/ace-s/ace_cem/' + pool_spot_counts_tsv_file_basename.replace('.tsv', '_ace-cem.tsv')
        pool.apply_async(
            run_ace_deconvolution,
            args=[
                pool_spot_counts_tsv_file,
                design_tsv_file,
                'cem',
                num_coverage,
                output_file
            ]
        )

    pool.close()
    pool.join()

    print("Finished successfully")
