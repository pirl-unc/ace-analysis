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
CONFIGURATIONS_DIR = "/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/02_benchmark_ace/02_benchmark_higher_positivity_rates/03_generated_configurations"
POOL_SPOT_COUNTS_DIR = "/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/02_benchmark_ace/02_benchmark_higher_positivity_rates/04_simulated_pool_spot_counts/300immunogenic_30nonimmunogenic_1dispersion_fnr0.00"
OUTPUT_DIR = "/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/02_benchmark_ace/02_benchmark_higher_positivity_rates/05_deconvolution_results/300immunogenic_30nonimmunogenic_1dispersion_fnr0.00"
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

    # Step 2. ACE-S configurations - ACE-EM
    if not os.path.exists(OUTPUT_DIR + '/ace-s/ace_em/'):
        os.makedirs(OUTPUT_DIR + '/ace-s/ace_em/')
    for pool_spot_counts_tsv_file in glob.glob(POOL_SPOT_COUNTS_DIR + '/ace-s/*pool_spot_counts.tsv'):
        pool_spot_counts_tsv_file_basename = os.path.basename(pool_spot_counts_tsv_file)
        num_peptides_per_pool = int(pool_spot_counts_tsv_file_basename.split('_')[3].replace('perpool',''))
        design_tsv_file = glob.glob(CONFIGURATIONS_DIR + '/ace-s/' + '_'.join(pool_spot_counts_tsv_file_basename.split('_')[0:7]) + '*.tsv')[0]
        output_file = OUTPUT_DIR + '/ace-s/ace_em/' + pool_spot_counts_tsv_file_basename.replace('.tsv', '_ace-em.tsv')
        pool.apply_async(
            run_ace_deconvolution,
            args=[
                pool_spot_counts_tsv_file,
                design_tsv_file,
                'em',
                3,
                output_file
            ]
        )

    # Step 3. ACE-S configurations - ACE-LASSO
    if not os.path.exists(OUTPUT_DIR + '/ace-s/ace_lasso/'):
        os.makedirs(OUTPUT_DIR + '/ace-s/ace_lasso/')
    for pool_spot_counts_tsv_file in glob.glob(POOL_SPOT_COUNTS_DIR + '/ace-s/*pool_spot_counts.tsv'):
        pool_spot_counts_tsv_file_basename = os.path.basename(pool_spot_counts_tsv_file)
        num_peptides_per_pool = int(pool_spot_counts_tsv_file_basename.split('_')[3].replace('perpool',''))
        design_tsv_file = glob.glob(CONFIGURATIONS_DIR + '/ace-s/' + '_'.join(pool_spot_counts_tsv_file_basename.split('_')[0:7]) + '*.tsv')[0]
        output_file = OUTPUT_DIR + '/ace-s/ace_lasso/' + pool_spot_counts_tsv_file_basename.replace('.tsv', '_ace-lasso.tsv')
        pool.apply_async(
            run_ace_deconvolution,
            args=[
                pool_spot_counts_tsv_file,
                design_tsv_file,
                'lasso',
                3,
                output_file
            ]
        )

    # Step 4. ACE-S configurations - ACE-EMPIRICAL
    if not os.path.exists(OUTPUT_DIR + '/ace-s/ace_empirical/'):
        os.makedirs(OUTPUT_DIR + '/ace-s/ace_empirical/')
    for pool_spot_counts_tsv_file in glob.glob(POOL_SPOT_COUNTS_DIR + '/ace-s/*pool_spot_counts.tsv'):
        pool_spot_counts_tsv_file_basename = os.path.basename(pool_spot_counts_tsv_file)
        num_peptides_per_pool = int(pool_spot_counts_tsv_file_basename.split('_')[3].replace('perpool',''))
        design_tsv_file = glob.glob(CONFIGURATIONS_DIR + '/ace-s/' + '_'.join(pool_spot_counts_tsv_file_basename.split('_')[0:7]) + '*.tsv')[0]
        output_file = OUTPUT_DIR + '/ace-s/ace_empirical/' + pool_spot_counts_tsv_file_basename.replace('.tsv', '_ace-empirical.tsv')
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

    # Step 5. ACE configurations - ACE-CEM
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

    # Step 6. ACE configurations - ACE-EM
    if not os.path.exists(OUTPUT_DIR + '/ace/ace_em/'):
        os.makedirs(OUTPUT_DIR + '/ace/ace_em/')
    for pool_spot_counts_tsv_file in glob.glob(POOL_SPOT_COUNTS_DIR + '/ace/*pool_spot_counts.tsv'):
        pool_spot_counts_tsv_file_basename = os.path.basename(pool_spot_counts_tsv_file)
        num_peptides_per_pool = int(pool_spot_counts_tsv_file_basename.split('_')[3].replace('perpool',''))
        design_tsv_file = glob.glob(CONFIGURATIONS_DIR + '/ace/' + '_'.join(pool_spot_counts_tsv_file_basename.split('_')[0:7]) + '*.tsv')[0]
        output_file = OUTPUT_DIR + '/ace/ace_em/' + pool_spot_counts_tsv_file_basename.replace('.tsv', '_ace-em.tsv')
        pool.apply_async(
            run_ace_deconvolution,
            args=[
                pool_spot_counts_tsv_file,
                design_tsv_file,
                'em',
                3,
                output_file
            ]
        )

    # Step 7. ACE configurations - ACE-LASSO
    if not os.path.exists(OUTPUT_DIR + '/ace/ace_lasso/'):
        os.makedirs(OUTPUT_DIR + '/ace/ace_lasso/')
    for pool_spot_counts_tsv_file in glob.glob(POOL_SPOT_COUNTS_DIR + '/ace/*pool_spot_counts.tsv'):
        pool_spot_counts_tsv_file_basename = os.path.basename(pool_spot_counts_tsv_file)
        num_peptides_per_pool = int(pool_spot_counts_tsv_file_basename.split('_')[3].replace('perpool',''))
        design_tsv_file = glob.glob(CONFIGURATIONS_DIR + '/ace/' + '_'.join(pool_spot_counts_tsv_file_basename.split('_')[0:7]) + '*.tsv')[0]
        output_file = OUTPUT_DIR + '/ace/ace_lasso/' + pool_spot_counts_tsv_file_basename.replace('.tsv', '_ace-lasso.tsv')
        pool.apply_async(
            run_ace_deconvolution,
            args=[
                pool_spot_counts_tsv_file,
                design_tsv_file,
                'lasso',
                3,
                output_file
            ]
        )

    # Step 8. ACE configurations - ACE-EMPIRICAL
    if not os.path.exists(OUTPUT_DIR + '/ace/ace_empirical/'):
        os.makedirs(OUTPUT_DIR + '/ace/ace_empirical/')
    for pool_spot_counts_tsv_file in glob.glob(POOL_SPOT_COUNTS_DIR + '/ace/*pool_spot_counts.tsv'):
        pool_spot_counts_tsv_file_basename = os.path.basename(pool_spot_counts_tsv_file)
        num_peptides_per_pool = int(pool_spot_counts_tsv_file_basename.split('_')[3].replace('perpool',''))
        design_tsv_file = glob.glob(CONFIGURATIONS_DIR + '/ace/' + '_'.join(pool_spot_counts_tsv_file_basename.split('_')[0:7]) + '*.tsv')[0]
        output_file = OUTPUT_DIR + '/ace/ace_empirical/' + pool_spot_counts_tsv_file_basename.replace('.tsv', '_ace-empirical.tsv')
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

    # Step 9. Repeated configurations - ACE-EMPIRICAL
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

    # Step 10. Repeated configurations - ACE-EM
    if not os.path.exists(OUTPUT_DIR + '/repeated/ace_em/'):
        os.makedirs(OUTPUT_DIR + '/repeated/ace_em/')
    for pool_spot_counts_tsv_file in glob.glob(POOL_SPOT_COUNTS_DIR + '/repeated/*pool_spot_counts.tsv'):
        pool_spot_counts_tsv_file_basename = os.path.basename(pool_spot_counts_tsv_file)
        num_peptides_per_pool = int(pool_spot_counts_tsv_file_basename.split('_')[3].replace('perpool',''))
        design_tsv_file = CONFIGURATIONS_DIR + '/repeated/' + '_'.join(pool_spot_counts_tsv_file_basename.split('_')[0:7]) + '.tsv'
        output_file = OUTPUT_DIR + '/repeated/ace_em/' + pool_spot_counts_tsv_file_basename.replace('.tsv', '_ace-em.tsv')
        pool.apply_async(
            run_ace_deconvolution,
            args=[
                pool_spot_counts_tsv_file,
                design_tsv_file,
                'em',
                3,
                output_file
            ]
        )

    # Step 11. Randomized configurations - ACE-CEM
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

    # Step 12. Randomized configurations - ACE-EM
    if not os.path.exists(OUTPUT_DIR + '/randomized/ace_em/'):
        os.makedirs(OUTPUT_DIR + '/randomized/ace_em/')
    tasks = []
    for pool_spot_counts_tsv_file in glob.glob(POOL_SPOT_COUNTS_DIR + '/randomized/*pool_spot_counts.tsv'):
        pool_spot_counts_tsv_file_basename = os.path.basename(pool_spot_counts_tsv_file)
        num_peptides_per_pool = int(pool_spot_counts_tsv_file_basename.split('_')[3].replace('perpool',''))
        design_tsv_file = glob.glob(CONFIGURATIONS_DIR + '/randomized/' + '_'.join(pool_spot_counts_tsv_file_basename.split('_')[0:7]) + '*.tsv')[0]
        output_file = OUTPUT_DIR + '/randomized/ace_em/' + pool_spot_counts_tsv_file_basename.replace('.tsv', '_ace-em.tsv')
        pool.apply_async(
            run_ace_deconvolution,
            args=[
                pool_spot_counts_tsv_file,
                design_tsv_file,
                'em',
                3,
                output_file
            ]
        )

    # Step 13. Randomized configurations - ACE-EMPIRICAL
    if not os.path.exists(OUTPUT_DIR + '/randomized/ace_empirical/'):
        os.makedirs(OUTPUT_DIR + '/randomized/ace_empirical/')
    tasks = []
    for pool_spot_counts_tsv_file in glob.glob(POOL_SPOT_COUNTS_DIR + '/randomized/*pool_spot_counts.tsv'):
        pool_spot_counts_tsv_file_basename = os.path.basename(pool_spot_counts_tsv_file)
        num_peptides_per_pool = int(pool_spot_counts_tsv_file_basename.split('_')[3].replace('perpool',''))
        design_tsv_file = glob.glob(CONFIGURATIONS_DIR + '/randomized/' + '_'.join(pool_spot_counts_tsv_file_basename.split('_')[0:7]) + '*.tsv')[0]
        output_file = OUTPUT_DIR + '/randomized/ace_empirical/' + pool_spot_counts_tsv_file_basename.replace('.tsv', '_ace-empirical.tsv')
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

    # Step 14. Randomized configurations - ACE-LASSO
    if not os.path.exists(OUTPUT_DIR + '/randomized/ace_lasso/'):
        os.makedirs(OUTPUT_DIR + '/randomized/ace_lasso/')
    tasks = []
    for pool_spot_counts_tsv_file in glob.glob(POOL_SPOT_COUNTS_DIR + '/randomized/*pool_spot_counts.tsv'):
        pool_spot_counts_tsv_file_basename = os.path.basename(pool_spot_counts_tsv_file)
        num_peptides_per_pool = int(pool_spot_counts_tsv_file_basename.split('_')[3].replace('perpool',''))
        design_tsv_file = glob.glob(CONFIGURATIONS_DIR + '/randomized/' + '_'.join(pool_spot_counts_tsv_file_basename.split('_')[0:7]) + '*.tsv')[0]
        output_file = OUTPUT_DIR + '/randomized/ace_lasso/' + pool_spot_counts_tsv_file_basename.replace('.tsv', '_ace-lasso.tsv')
        pool.apply_async(
            run_ace_deconvolution,
            args=[
                pool_spot_counts_tsv_file,
                design_tsv_file,
                'lasso',
                3,
                output_file
            ]
        )

    pool.close()
    pool.join()
    print("Finished successfully")
