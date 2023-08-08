"""
The purpose of this python3 script is to implement the Solver dataclass
"""


from dataclasses import dataclass, field
from .solver import Solver
from acelib.main import run_ace_generate
from acelib.block_assignment import BlockAssignment
from acelib.types import *


@dataclass(frozen=True)
class AceSolver(Solver):
    cluster_peptides: bool
    mode: str
    trained_model_file: str
    sim_threshold: float = 0.8
    sim_fxn: str = 'euclidean'
    golfy_max_iters: int = 2000
    golfy_init_mode: str = 'greedy'
    golfy_allow_extra_pools: bool = False
    max_peptides_per_block: int = 100
    max_peptides_per_pool: int = 10
    num_processes: int = 1
    shuffle_iters: int = 1000
    assign_well_ids: bool = True
    num_plate_wells = 96
    verbose: bool = False

    def generate_assignment(
            self,
            peptides: Peptides,
            num_peptides_per_pool: int,
            num_coverage: int,
            random_seed: int
    ) -> Tuple[BlockAssignment, PeptidePairs]:
        """
        Generates an ELISpot configuration.

        Parameters
        ----------
        peptides                    :   Peptides (list of tuples (peptide ID, peptide sequence)).
        num_peptides_per_pool       :   Number of peptides per pool
        num_coverage                :   Number of coverage.
        random_seed                 :   Random seed.    
        
        Returns
        -------
        block_assignment            :   BlockAssignment object.
        preferred_peptide_pairs     :   PeptidePairs.
        """
        block_assignment, block_design = run_ace_generate(
            peptides=peptides,
            num_peptides_per_pool=num_peptides_per_pool,
            num_coverage=num_coverage,
            cluster_peptides=self.cluster_peptides,
            trained_model_file=self.trained_model_file,
            mode=self.mode,
            sequence_similarity_function=self.sim_fxn,
            sequence_similarity_threshold=self.sim_threshold,
            golfy_random_seed=random_seed,
            golfy_strategy=self.golfy_init_mode,
            golfy_max_iters=self.golfy_max_iters,
            golfy_allow_extra_pools=self.golfy_allow_extra_pools,
            cpsat_solver_num_processes=self.num_processes,
            cpsat_solver_shuffle_iters=self.shuffle_iters,
            cpsat_solver_max_peptides_per_block=self.max_peptides_per_block,
            cpsat_solver_max_peptides_per_pool=self.max_peptides_per_pool,
            assign_well_ids=self.assign_well_ids,
            num_plate_wells=self.num_plate_wells,
            verbose=False
        )
        return block_assignment, block_design.preferred_peptide_pairs
