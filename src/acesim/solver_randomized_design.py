"""
The purpose of this python3 script is to implement the Bogey dataclass
"""


import random
import pandas as pd
from dataclasses import dataclass, field
from acelib.block_assignment import BlockAssignment
from acelib.utilities import convert_peptides_to_dataframe
from acelib.types import *
from .solver import Solver


@dataclass(frozen=True)
class RandomizedDesignSolver(Solver):
    random_seed: int

    def generate_assignment(
            self,
            peptides: Peptides,
            num_peptides_per_pool: int,
            num_coverage: int
    ) -> Tuple[BlockAssignment, PeptidePairs]:
        """
        Randomly generates an ELISpot configuration.

        Parameters
        ----------
        peptides                    :   Peptides.
        num_peptides_per_pool       :   Number of peptides per pool
        num_coverage                :   Number of coverage.

        Returns
        -------
        block_assignment            :   BlockAssignment object.
        preferred_peptide_pairs     :   PeptidePairs.
        """
        random.seed(self.random_seed)
        df_peptides = convert_peptides_to_dataframe(peptides=peptides)

        # Step 1. Initialize pools
        pools = {}
        num_peptides = len(df_peptides['peptide_id'].unique())
        num_pools_per_coverage = int(num_peptides / num_peptides_per_pool)

        for coverage in range(1, num_coverage + 1):
            pools[coverage] = {}
            for j in range(1, num_pools_per_coverage + 1):
                pool = num_pools_per_coverage * (coverage - 1) + j
                pools[coverage][pool] = []

        # Step 2. Randomly assign peptides to pools
        peptide_ids = list(df_peptides['peptide_id'].unique())
        for coverage in range(1, num_coverage + 1):
            peptide_ids_ = peptide_ids.copy()
            random.shuffle(peptide_ids_)
            for peptide_id in peptide_ids_:
                for pool in pools[coverage].keys():
                    if len(pools[coverage][pool]) < num_peptides_per_pool:
                        pools[coverage][pool].append(peptide_id)
                        break

        # Step 3. Create a BlockAssignment object
        block_assignment = BlockAssignment()
        for coverage, value in pools.items():
            for pool, value2 in value.items():
                for peptide_id in value2:
                    df_curr_peptide = df_peptides.loc[df_peptides['peptide_id'] == peptide_id,:]
                    block_assignment.add_peptide(
                        coverage=coverage,
                        pool=pool,
                        peptide_id=peptide_id,
                        peptide_sequence=df_curr_peptide['peptide_sequence'].values[0]
                    )
        return block_assignment, []
