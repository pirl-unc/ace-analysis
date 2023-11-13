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


NUM_PROCESSES = 32
CONFIGURATIONS_DIR = "/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/04_real_world_datasets/tan_2021/03_generated_configurations"
POOL_SPOT_COUNTS_DIR = "/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/04_real_world_datasets/tan_2021/04_simulated_pool_spot_counts/300immunogenic_30nonimmunogenic_1dispersion_fnr0.00"
OUTPUT_DIR = "/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/04_real_world_datasets/tan_2021/05_deconvolution_results/300immunogenic_30nonimmunogenic_1dispersion_fnr0.00"
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
    pool = mp.Pool(processes=NUM_PROCESSES)

    # Step 1. ACE-S configurations - ACE-CEM
    if not os.path.exists(OUTPUT_DIR + '/ace-s/ace_cem/'):
        os.makedirs(OUTPUT_DIR + '/ace-s/ace_cem/')
    for pool_spot_counts_tsv_file in glob.glob(POOL_SPOT_COUNTS_DIR + '/ace-s/*pool_spot_counts.tsv'):
        pool_spot_counts_tsv_file_basename = os.path.basename(pool_spot_counts_tsv_file)
        num_peptides_per_pool = int(pool_spot_counts_tsv_file_basename.split('_')[3].replace('perpool',''))
        design_tsv_file = glob.glob(CONFIGURATIONS_DIR + '/ace-s/' + '_'.join(pool_spot_counts_tsv_file_basename.split('_')[0:7]) + '*.tsv')[0]
        output_file = OUTPUT_DIR + '/ace-s/ace_cem/' + pool_spot_counts_tsv_file_basename.replace('.tsv', '_ace-cem.tsv')
        pool.apply_async(
            run_ace_deconvolution,
            args=[
                pool_spot_counts_tsv_file,
                design_tsv_file,
                'cem',
                3,
                output_file
            ]
        )

    # Step 2. ACE configurations - ACE-CEM
    if not os.path.exists(OUTPUT_DIR + '/ace/ace_cem/'):
        os.makedirs(OUTPUT_DIR + '/ace/ace_cem/')
    for pool_spot_counts_tsv_file in glob.glob(POOL_SPOT_COUNTS_DIR + '/ace/*pool_spot_counts.tsv'):
        pool_spot_counts_tsv_file_basename = os.path.basename(pool_spot_counts_tsv_file)
        num_peptides_per_pool = int(pool_spot_counts_tsv_file_basename.split('_')[3].replace('perpool',''))
        design_tsv_file = glob.glob(CONFIGURATIONS_DIR + '/ace/' + '_'.join(pool_spot_counts_tsv_file_basename.split('_')[0:7]) + '*.tsv')[0]
        output_file = OUTPUT_DIR + '/ace/ace_cem/' + pool_spot_counts_tsv_file_basename.replace('.tsv', '_ace-cem.tsv')
        pool.apply_async(
            run_ace_deconvolution,
            args=[
                pool_spot_counts_tsv_file,
                design_tsv_file,
                'cem',
                3,
                output_file
            ]
        )

    # Step 3. Repeated configurations - ACE-EMPIRICAL
    if not os.path.exists(OUTPUT_DIR + '/repeated/ace_empirical/'):
        os.makedirs(OUTPUT_DIR + '/repeated/ace_empirical/')
    for pool_spot_counts_tsv_file in glob.glob(POOL_SPOT_COUNTS_DIR + '/repeated/*pool_spot_counts.tsv'):
        pool_spot_counts_tsv_file_basename = os.path.basename(pool_spot_counts_tsv_file)
        num_peptides_per_pool = int(pool_spot_counts_tsv_file_basename.split('_')[3].replace('perpool',''))
        design_tsv_file = CONFIGURATIONS_DIR + '/repeated/' + '_'.join(pool_spot_counts_tsv_file_basename.split('_')[0:7]) + '.tsv'
        output_file = OUTPUT_DIR + '/repeated/ace_empirical/' + pool_spot_counts_tsv_file_basename.replace('.tsv', '_ace-empirical.tsv')
        pool.apply_async(
            run_ace_deconvolution,
            args=[
                pool_spot_counts_tsv_file,
                design_tsv_file,
                'empirical',
                3,
                output_file
            ]
        )

    # Step 4. Randomized configurations - ACE-CEM
    if not os.path.exists(OUTPUT_DIR + '/randomized/ace_cem/'):
        os.makedirs(OUTPUT_DIR + '/randomized/ace_cem/')
    tasks = []
    for pool_spot_counts_tsv_file in glob.glob(POOL_SPOT_COUNTS_DIR + '/randomized/*pool_spot_counts.tsv'):
        pool_spot_counts_tsv_file_basename = os.path.basename(pool_spot_counts_tsv_file)
        num_peptides_per_pool = int(pool_spot_counts_tsv_file_basename.split('_')[3].replace('perpool',''))
        design_tsv_file = glob.glob(CONFIGURATIONS_DIR + '/randomized/' + '_'.join(pool_spot_counts_tsv_file_basename.split('_')[0:7]) + '*.tsv')[0]
        output_file = OUTPUT_DIR + '/randomized/ace_cem/' + pool_spot_counts_tsv_file_basename.replace('.tsv', '_ace-cem.tsv')
        pool.apply_async(
            run_ace_deconvolution,
            args=[
                pool_spot_counts_tsv_file,
                design_tsv_file,
                'cem',
                3,
                output_file
            ]
        )

    pool.close()
    pool.join()
