"""
The purpose of this python3 script is to implement the Solver dataclass
"""


from dataclasses import dataclass, field
from .configurator import Configurator
from acelib.main import run_ace_generate
from acelib.block_assignment import BlockAssignment
from acelib.types import *


@dataclass(frozen=True)
class AceConfigurator(Configurator):

    def generate_assignment(
            self,
            peptides: Peptides,
            num_peptides_per_pool: int,
            num_coverage: int,
            trained_model_file: str,
            cluster_peptides: bool,
            mode: str,
            sequence_similarity_function: str,
            sequence_similarity_threshold: float,
            golfy_random_seed: int,
            golfy_strategy: str,
            golfy_max_iters: int,
            golfy_allow_extra_pools: bool,
            plate_size: int
    ) -> Tuple[BlockAssignment, PreferredPeptidePairs]:
        block_assignment, block_design = run_ace_generate(
            peptides=peptides,
            num_peptides_per_pool=num_peptides_per_pool,
            num_coverage=num_coverage,
            trained_model_file=trained_model_file,
            cluster_peptides=cluster_peptides,
            mode=mode,
            sequence_similarity_function=sequence_similarity_function,
            sequence_similarity_threshold=sequence_similarity_threshold,
            golfy_random_seed=golfy_random_seed,
            golfy_strategy=golfy_strategy,
            golfy_max_iters=golfy_max_iters,
            golfy_allow_extra_pools=golfy_allow_extra_pools,
            plate_size=plate_size,
            verbose=False
        )
        return block_assignment, block_design.preferred_peptide_pairs
