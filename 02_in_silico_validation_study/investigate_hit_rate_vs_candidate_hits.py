"""
The purpose of this python3 is to perform in silico validation experiment on
investigation of hit count vs number of candidate hits.
"""


import glob
import multiprocessing as mp
import pandas as pd
import random
import os
from typing import Dict
from acelib.constants import DeconvolutionResults
from acelib.main import run_ace_identify
from functools import partial


MAX_POSITIVE_PEPTIDE_HIT_RATE = 0.5 # proportion of all peptides that are (ground truth) positive hits
NUM_REPEATS = 10000
NUM_THREADS = 48
DATA_DIR = '/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data'
OUTPUT_DIR = DATA_DIR + '/processed/02_in_silico_validation_study'


def run_randomized_hit_experiment(
        df_configuration: pd.DataFrame,
        hit_count: int,
        num_initial_pools: int,
        experiment_id: int
) -> Dict:
    """
    Simulates hit peptides

    Parameters
    ----------
    df_configuration    :   DataFrame with the following columns:
                            'coverage_id'
                            'pool_id'
                            'peptide_id'
    hit_count           :   Ground truth hit count.
    num_initial_pools   :   Number of initial pools.
    experiment_id       :   Experiment ID.

    Returns
    -------
    results_dict        :   Dictionary with the following keys:
                            'ground_truth_hit_peptide_ids'
                            'ground_truth_hit_peptide_ids_count'
                            'predicted_hit_peptide_ids'
                            'predicted_hit_peptide_ids_count'
                            'candidate_hit_peptide_ids'
                            'candidate_hit_peptide_ids_count'
                            'num_initial_pools'
                            'num_total_pools'
    """
    data = {}

    # Step 1. Randomly select peptide(s)
    peptide_ids = list(df_configuration['peptide_id'].unique())
    ground_truth_hit_peptide_ids = random.sample(peptide_ids, hit_count)

    # Step 2. Identify hit pool IDs associated with randomly selected peptides
    hit_pool_ids = list(df_configuration.loc[df_configuration['peptide_id'].isin(ground_truth_hit_peptide_ids), 'pool_id'].unique())

    # Step 3. Identify hit peptide IDs
    df_hits_max = run_ace_identify(
        hit_pool_ids=hit_pool_ids,
        df_configuration=df_configuration
    )

    # Step 4. Identify predicted hit peptide IDs
    predicted_hit_peptide_ids = df_hits_max.loc[df_hits_max['deconvolution_result'] == DeconvolutionResults.HIT, 'peptide_id'].values.tolist()

    # Step 5. Check that each predicted peptide ID is in the list of ground truth hit peptide IDs
    for peptide_id in predicted_hit_peptide_ids:
        if peptide_id not in ground_truth_hit_peptide_ids:
            print('ERROR: %s is not in the original set of ground truth hit peptide IDs.' % peptide_id)
            exit(1)

    # Step 6. Identify candidate hit peptide IDs
    candidate_hit_peptide_ids = df_hits_max.loc[df_hits_max['deconvolution_result'] == DeconvolutionResults.CANDIDATE_HIT, 'peptide_id'].values.tolist()

    # Step 7. Store experiment results
    data['experiment_id'] = experiment_id
    data['ground_truth_hit_peptide_ids'] = ';'.join(ground_truth_hit_peptide_ids)
    data['ground_truth_hit_peptide_ids_count'] = len(ground_truth_hit_peptide_ids)
    data['predicted_hit_peptide_ids'] = ';'.join(predicted_hit_peptide_ids)
    data['predicted_hit_peptide_ids_count'] = len(predicted_hit_peptide_ids)
    data['candidate_hit_peptide_ids'] = ';'.join(candidate_hit_peptide_ids)
    data['candidate_hit_peptide_ids_count'] = len(candidate_hit_peptide_ids)
    data['num_initial_pools'] = num_initial_pools
    data['num_total_pools'] = num_initial_pools + len(candidate_hit_peptide_ids)
    return data


if __name__ == "__main__":
    for configuration_csv_file in glob.glob(DATA_DIR + '/processed/01_resource_configuration_generation/*.csv'):
        csv_file_basename = os.path.basename(configuration_csv_file)
        df_configuration = pd.read_csv(configuration_csv_file)

        # Figure out the initial number of pools
        max_pool_idx = -1
        for pool_id in df_configuration['pool_id'].unique():
            pool_idx = int(pool_id.split('_')[1])
            if pool_idx > max_pool_idx:
                max_pool_idx = pool_idx

        data = {
            'experiment_id': [],
            'ground_truth_hit_peptide_ids': [],
            'ground_truth_hit_peptide_ids_count': [],
            'predicted_hit_peptide_ids': [],
            'predicted_hit_peptide_ids_count': [],
            'candidate_hit_peptide_ids': [],
            'candidate_hit_peptide_ids_count': [],
            'num_initial_pools': [],
            'num_total_pools': []
        }
        num_peptides = len(df_configuration['peptide_id'].unique())
        curr_hit_count = 1
        curr_hit_rate = curr_hit_count / num_peptides
        curr_experiment_id = 1
        while curr_hit_rate <= MAX_POSITIVE_PEPTIDE_HIT_RATE:
            pool = mp.Pool(processes=NUM_THREADS)
            func = partial(run_randomized_hit_experiment, df_configuration, curr_hit_count, max_pool_idx)
            results = pool.map(func, [i for i in range(1, NUM_REPEATS + 1)])
            pool.close()
            for result in results:
                data['experiment_id'].append(result['experiment_id'])
                data['ground_truth_hit_peptide_ids'].append(result['ground_truth_hit_peptide_ids'])
                data['ground_truth_hit_peptide_ids_count'].append(result['ground_truth_hit_peptide_ids_count'])
                data['predicted_hit_peptide_ids'].append(result['predicted_hit_peptide_ids'])
                data['predicted_hit_peptide_ids_count'].append(result['predicted_hit_peptide_ids_count'])
                data['candidate_hit_peptide_ids'].append(result['candidate_hit_peptide_ids'])
                data['candidate_hit_peptide_ids_count'].append(result['candidate_hit_peptide_ids_count'])
                data['num_initial_pools'].append(result['num_initial_pools'])
                data['num_total_pools'].append(result['num_total_pools'])
            curr_hit_count += 1
            curr_hit_rate = curr_hit_count / num_peptides

        df_results = pd.DataFrame(data)
        df_results.to_csv(OUTPUT_DIR + '/' + csv_file_basename + '_hit_counts_vs_candidate_hits.tsv', sep='\t', index=False)

