"""
The purpose of this python3 script is to implement the Solver dataclass
"""


import pandas as pd
from dataclasses import dataclass, field
from .solver import Solver


@dataclass(frozen=True)
class AceSolver(Solver):

    def generate(
            self,
            df_peptides: pd.DataFrame
    ):
        # todo run ace
        pass

