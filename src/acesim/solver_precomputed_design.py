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
class PrecomputedSolver(Solver):
    block_assignment: BlockAssignment

    def generate_assignment(
            self,
            peptides: Peptides,
            num_peptides_per_pool: int,
            num_coverage: int,
            random_seed: int
    ) -> Tuple[BlockAssignment, PeptidePairs]:
        return self.block_assignment, []
