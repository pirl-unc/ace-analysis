"""
The purpose of this python3 script is to implement the AceDeconvolver dataclass
"""


import pandas as pd
from dataclasses import dataclass, field
from acelib.main import run_ace_deconvolve
from acelib.block_assignment import BlockAssignment
from acelib.deconvolution import DeconvolutionResult
from .deconvolver import Deconvolver


@dataclass(frozen=True)
class PrecomputedDeconvolver(Deconvolver):
    deconvolution_result: DeconvolutionResult

    def deconvolve(self) -> DeconvolutionResult:
        return self.deconvolution_result
