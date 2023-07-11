"""
The purpose of this python3 script is to implement the Bogey dataclass
"""


import random
import pandas as pd
from dataclasses import dataclass, field
from .solver import Solver


@dataclass(frozen=True)
class BogeySolver(Solver):
    random_seed: int

    def generate_configuration(
            self,
            df_peptides: pd.DataFrame,
            num_peptides_per_pool: int,
            num_coverage: int
    ) -> pd.DataFrame:
        """
        Randomly generates an ELISpot configuration.

        Parameters
        ----------
        df_peptides                 :   pd.DataFrame with the following columns:
                                        'peptide_id'
                                        'peptide_sequence'
        num_peptides_per_pool       :   Number of peptides per pool
        num_coverage                :   Number of coverage.

        Returns
        -------
        df_configuration    :   pd.DataFrame with the following columns:
                                'coverage_id'
                                'pool_id'
                                'peptide_id'
                                'peptide_sequence'
        """
        random.seed(self.random_seed)

        # Step 1. Initialize pools
        pools = {}
        num_peptides = len(df_peptides['peptide_id'].unique())
        num_pools_per_coverage = int(num_peptides / num_peptides_per_pool)

        for i in range(1, num_coverage + 1):
            pools['coverage_%i' % i] = {}
            for j in range(1, num_pools_per_coverage + 1):
                pool_id = 'pool_%i' % (num_pools_per_coverage * (i - 1) + j)
                pools['coverage_%i' % i][pool_id] = []

        # Step 2. Randomly assign peptides to pools
        peptide_ids = list(df_peptides['peptide_id'].unique())
        for i in range(1, num_coverage + 1):
            peptide_ids_ = peptide_ids.copy()
            random.shuffle(peptide_ids_)
            for peptide_id in peptide_ids_:
                for pool_id in pools['coverage_%i' % i].keys():
                    if len(pools['coverage_%i' % i][pool_id]) < num_peptides_per_pool:
                        pools['coverage_%i' % i][pool_id].append(peptide_id)
                        break

        # Step 3. Convert assignments to a DataFrame
        data = {
            'coverage_id': [],
            'pool_id': [],
            'peptide_id': [],
            'peptide_sequence': []
        }
        for coverage_id, value in pools.items():
            for pool_id, value2 in value.items():
                for peptide_id in value2:
                    df_curr_peptide = df_peptides.loc[df_peptides['peptide_id'] == peptide_id,:]
                    data['coverage_id'].append(coverage_id)
                    data['pool_id'].append(pool_id)
                    data['peptide_id'].append(peptide_id)
                    data['peptide_sequence'].append(df_curr_peptide['peptide_sequence'].values[0])
        return pd.DataFrame(data)
