"""
The purpose of this python3 script is to implement the AceDeconvolver dataclass
"""


import pandas as pd
from typing import Literal
from dataclasses import dataclass, field
from acelib.main import run_ace_deconvolve
from acelib.block_assignment import BlockAssignment
from acelib.deconvolution import DeconvolutionResult
from acelib.constants import DeconvolutionMethods
from .deconvolver import Deconvolver


@dataclass(frozen=True)
class AceDeconvolver(Deconvolver):

    def deconvolve(
            self,
            df_readout: pd.DataFrame,
            block_assignment: BlockAssignment,
            method: str,
            min_coverage: int,
            min_pool_spot_count: float
    ) -> DeconvolutionResult:
        # Step 1. Perform empirical deconvolution
        deconvolution_result = run_ace_deconvolve(
            df_readout=df_readout,
            block_assignment=block_assignment,
            method=method,
            min_coverage=min_coverage,
            min_pool_spot_count=min_pool_spot_count,
            verbose=False
        )
        return deconvolution_result
