"""
The purpose of this python3 script is to implement the Solver dataclass
"""


import pandas as pd
import torch
from dataclasses import dataclass, field
from .solver import Solver
from acelib.block_assignment import BlockAssignment
from acelib.block_design import BlockDesign
from acelib.main import run_ace_golfy, run_ace_sat_solver
from acelib.sequence_features import AceNeuralEngine
from acelib.types import *
from transformers import AutoTokenizer, AutoModelForMaskedLM


@dataclass(frozen=True)
class AceSolver(Solver):
    cluster_peptides: bool
    random_seed: int
    mode: str
    trained_model_file: str
    sim_threshold: float = 0.7
    sim_fxn: str = 'euclidean'
    golfy_max_iters: int = 2000
    golfy_init_mode: str = 'greedy'
    golfy_allow_extra_pools: bool = False
    max_peptides_per_block = 100
    max_peptides_per_pool = 10
    num_processes = 1

    def generate_assignment(
            self,
            peptides: Peptides,
            num_peptides_per_pool: int,
            num_coverage: int
    ) -> Tuple[BlockAssignment, PeptidePairs]:
        """
        Generates an ELISpot configuration.

        Parameters
        ----------
        peptides                    :   Peptides (list of tuples (peptide ID, peptide sequence)).
        num_peptides_per_pool       :   Number of peptides per pool
        num_coverage                :   Number of coverage.
        
        Returns
        -------
        block_assignment            :   BlockAssignment object.
        preferred_peptide_pairs     :   PeptidePairs.
        """
        # Step 1. Identify pairs of similar peptides
        if self.cluster_peptides:
            # Load model
            ESM2_TOKENIZER = AutoTokenizer.from_pretrained("facebook/esm2_t6_8M_UR50D")
            ESM2_MODEL = AutoModelForMaskedLM.from_pretrained("facebook/esm2_t6_8M_UR50D", return_dict=True, output_hidden_states=True)
            device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
            ace_eng = AceNeuralEngine(ESM2_MODEL, ESM2_TOKENIZER, device)
            ace_eng.load_weights(self.trained_model_file)

            # Predict similar peptides
            preferred_peptide_pairs = ace_eng.find_paired_peptides(
                peptide_ids=[p[0] for p in peptides],
                peptide_sequences=[p[1] for p in peptides],
                sim_fxn=self.sim_fxn,
                threshold=self.sim_threshold
            )
        else:
            preferred_peptide_pairs = []
        preferred_peptide_pairs = [(p1, p2) for p1, p2, score in preferred_peptide_pairs]

        # Step 2. Generate a block design
        block_design = BlockDesign(
            peptides=peptides,
            num_peptides_per_pool=num_peptides_per_pool,
            num_coverage=num_coverage,
            max_peptides_per_block=self.max_peptides_per_block,
            disallowed_peptide_pairs=[],
            preferred_peptide_pairs=preferred_peptide_pairs
        )

        # Step 3. Generate a block assignment
        if self.mode == 'golfy':
            block_assignment = run_ace_golfy(
                block_design=block_design,
                random_seed=self.random_seed,
                max_iters=self.golfy_max_iters,
                init_mode=self.golfy_init_mode,
                allow_extra_pools=self.golfy_allow_extra_pools,
                verbose=False
            )
        elif self.mode == 'sat_solver':
            block_assignment = run_ace_sat_solver(
                block_design=block_design,
                max_peptides_per_pool=self.max_peptides_per_pool,
                num_processes=self.num_processes,
                verbose=True
            )
        else:
            print("Unknown mode: %s" % self.mode)
            exit(1)
        return block_assignment, preferred_peptide_pairs
