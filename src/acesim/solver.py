"""
The purpose of this python3 script is to implement the Solver dataclass
"""


import pandas as pd
from dataclasses import dataclass


@dataclass(frozen=True)
class Solver:
    name: str

    def generate_configuration(
            self,
            df_peptides: pd.DataFrame,
            num_peptides_per_pool: int,
            num_coverage: int,
            **kwargs
    ) -> pd.DataFrame:
        raise Exception("Subclass must implement 'generate_configuration' method")

