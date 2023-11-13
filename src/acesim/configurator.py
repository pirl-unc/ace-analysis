"""
The purpose of this python3 script is to implement the Configurator dataclass
"""


import pandas as pd
from dataclasses import dataclass
from acelib.block_assignment import BlockAssignment
from acelib.types import *


@dataclass(frozen=True)
class Configurator:
    name: str

    def generate_assignment(self, **kwargs) -> Tuple[BlockAssignment, PreferredPeptidePairs]:
        raise Exception("Subclass must implement 'generate_assignment' method")

