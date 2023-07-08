"""
The purpose of this python3 script is to implement the Experiment dataclass
"""


from dataclasses import dataclass, field
from .solver import Solver


@dataclass(frozen=True)
class Experiment:
    solver: Solver

    def run(self):
        pass
