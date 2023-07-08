"""
The purpose of this python3 script is to implement the Solver dataclass
"""


import pandas as pd
from dataclasses import dataclass, field


@dataclass(frozen=True)
class Solver:

    def generate(
            self,
            df_peptides: pd.DataFrame,
            num_peptides_per_pool: int
    ):
        raise Exception("Subclass must implement 'generate' method")

