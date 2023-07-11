"""
The purpose of this python3 script is to implement the Solver dataclass
"""


import pandas as pd
import torch
from dataclasses import dataclass, field
from .solver import Solver
from acelib.main import run_ace_golfy, run_ace_sat_solver
from acelib.sequence_features import AceNeuralEngine
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
    num_processes_per_batch = 100
    num_processes = 1

    def generate_configuration(
            self,
            df_peptides: pd.DataFrame,
            num_peptides_per_pool: int,
            num_coverage: int
    ) -> pd.DataFrame:
        """
        Generates an ELISpot configuration.

        Parameters
        ----------
        df_peptides                 :   pd.DataFrame with the following columns:
                                        'peptide_id'
                                        'peptide_sequence'
        num_peptides_per_pool       :   Number of peptides per pool
        num_coverage                :   Number of coverage.
        
        Returns
        -------
        df_configuration            :   pd.DataFrame with the following columns:
                                        'coverage_id'
                                        'pool_id'
                                        'peptide_id'
                                        'peptide_sequence'
        """
        if self.mode == 'golfy':
            if self.cluster_peptides:
                # Load model
                ESM2_TOKENIZER = AutoTokenizer.from_pretrained("facebook/esm2_t6_8M_UR50D")
                ESM2_MODEL = AutoModelForMaskedLM.from_pretrained("facebook/esm2_t6_8M_UR50D", return_dict=True, output_hidden_states=True)
                device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
                ace_eng = AceNeuralEngine(ESM2_MODEL, ESM2_TOKENIZER, device)
                ace_eng.load_weights(self.trained_model_file)

                # Predict similar peptides
                preferred_peptide_pairs = ace_eng.find_paired_peptides(
                    peptide_ids=df_peptides['peptide_id'].values.tolist(),
                    peptide_sequences=df_peptides['peptide_sequence'].values.tolist(),
                    sim_fxn=self.sim_fxn,
                    threshold=self.sim_threshold
                )

                # Run golfy
                is_valid, df_configuration = run_ace_golfy(
                    df_peptides=df_peptides,
                    num_peptides_per_pool=num_peptides_per_pool,
                    num_coverage=num_coverage,
                    random_seed=self.random_seed,
                    max_iters=self.golfy_max_iters,
                    init_mode=self.golfy_init_mode,
                    preferred_peptide_pairs=preferred_peptide_pairs
                )
            else:
                # Run golfy
                is_valid, df_configuration = run_ace_golfy(
                    df_peptides=df_peptides,
                    num_peptides_per_pool=num_peptides_per_pool,
                    num_coverage=num_coverage,
                    random_seed=self.random_seed,
                    max_iters=self.golfy_max_iters,
                    init_mode=self.golfy_init_mode
                )
            return df_configuration
        elif self.mode == 'sat_solver':
            if self.cluster_peptides:
                # Load model
                ESM2_TOKENIZER = AutoTokenizer.from_pretrained("facebook/esm2_t6_8M_UR50D")
                ESM2_MODEL = AutoModelForMaskedLM.from_pretrained("facebook/esm2_t6_8M_UR50D", return_dict=True, output_hidden_states=True)
                device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
                ace_eng = AceNeuralEngine(ESM2_MODEL, ESM2_TOKENIZER, device)
                ace_eng.load_weights(self.trained_model_file)

                # Predict similar peptides
                preferred_peptide_pairs = ace_eng.find_paired_peptides(
                    peptide_ids=df_peptides['peptide_id'].values.tolist(),
                    peptide_sequences=df_peptides['peptide_sequence'].values.tolist(),
                    sim_fxn=self.sim_fxn,
                    threshold=self.sim_threshold
                )

                # Run SAT solver
                df_configuration = run_ace_sat_solver(
                    df_peptides=df_peptides,
                    num_peptides_per_pool=num_peptides_per_pool,
                    num_coverage=num_coverage,
                    num_peptides_per_batch=self.num_peptides_per_batch,
                    random_seed=self.random_seed,
                    num_processes=self.num_processes,
                    preferred_peptide_pairs=preferred_peptide_pairs
                )
            else:
                # Run SAT solver
                df_configuration = run_ace_sat_solver(
                    df_peptides=df_peptides,
                    num_peptides_per_pool=num_peptides_per_pool,
                    num_coverage=num_coverage,
                    num_peptides_per_batch=self.num_peptides_per_batch,
                    random_seed=self.random_seed,
                    num_processes=self.num_processes
                )
            return df_configuration
        else:
            print("Unknown mode: %s" % self.mode)
            exit(1)
    
