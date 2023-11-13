"""
The purpose of this python3 script is to implement the Deconvolver dataclass
"""


import pandas as pd
from dataclasses import dataclass
from acelib.block_assignment import BlockAssignment
from acelib.deconvolution import DeconvolutionResult


@dataclass(frozen=True)
class Deconvolver:
    name: str

    def deconvolve(
            self,
            df_readout: pd.DataFrame,
            block_assignment: BlockAssignment,
            num_coverage: int,
            **kwargs
    ) -> DeconvolutionResult:
        raise Exception("Subclass must implement 'deconvolve' method")
