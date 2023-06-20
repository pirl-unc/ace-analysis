"""
The purpose of this python3 file is to implement 
baseline methods (i.e. benchmark approaches) for in silico validation.
"""


import numpy as np
import pandas as pd
import random
from typing import Tuple, List


def generate_random_elispot_configuration(
    peptide_ids: List[str],
    peptides_per_pool: int,
    num_coverage: int) -> pd.DataFrame:
    """
    Generates a random ELIspot configuration.

    Parameters
    ----------
    peptide_ids         :   List of peptide IDs.
    peptides_per_pool   :   Number of peptides per pool.
    num_coverage        :   Coverage number (i.e. repeats per peptide).

    Returns
    -------
    DataFrame with the following columns:
    'peptide_id'
    'pool_id'
    'coverage_id'
    """
    data = {
        'peptide_id': [],
        'pool_id': [],
        'coverage_id': []
    }
    for curr_coverage in range(1, num_coverage + 1):
        peptide_ids_temp = list(np.copy(peptide_ids))
        num_pools = int(len(peptide_ids) / peptides_per_pool)
        pool_ids = ['pool_' + str(i + ((curr_coverage - 1) * num_pools)) for i in range(1, num_pools + 1)]
        for curr_pool_id in pool_ids:
            random_peptide_ids = random.sample(peptide_ids_temp, k=peptides_per_pool)
            for curr_peptide_id in random_peptide_ids:
                data['peptide_id'].append(curr_peptide_id)
                data['pool_id'].append(curr_pool_id)
                data['coverage_id'].append(curr_coverage)
                peptide_ids_temp.remove(curr_peptide_id)
    df = pd.DataFrame(data)
    return df

