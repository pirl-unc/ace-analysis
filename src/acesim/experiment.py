"""
The purpose of this python3 script is to implement the Experiment dataclass
"""


import copy
import Levenshtein
import math
import multiprocessing as mp
import numpy as np
import pandas as pd
import random
from sklearn.metrics import roc_auc_score
from dataclasses import dataclass, field
from functools import partial
from typing import Dict, List, Literal, Tuple
from acelib.constants import DeconvolutionLabels, DeconvolveModes
from acelib.main import run_ace_deconvolve
from acelib.utilities import convert_dataframe_to_peptides
from .solver import Solver
from .solver_precomputed_design import PrecomputedSolver


@dataclass
class Experiment:
    """
    Class to simulate an end-to-end ELISpot experiment (generation of config,
    simulation of spot counts, and identification of positive peptides).
    Class accepts a list of Solver objects, which are used to generate a
    configuration for a given sample of peptides. For each sample of peptides,
    the class simulates the number of spots per peptide and maps this to each
    well per configuration from each solver. Then, the putative positive
    peptides are identified and test statistics are calculated for the
    deconvolution results.

    Simulation Procedure:
    1. Initialize the object with the given assay parameters.
    2. Load the reference database of peptides to sample from
    3. Generate a configuration of peptides to be used in the experiment
    4. Simulate the spot counts for each well in the experiment
    5. Identify the hit peptides from the simulated spot counts
    6. Use the candidate + hit peptides list to run the deconvolution
    7. Return statistics on the deconvolution results
    """
    # Experiment ID
    experiment_id: int
    
    # Define the assay parameters
    num_peptides: int
    num_positives: int
    num_peptides_per_pool: int
    coverage: int
    
    # Define the solver parameters
    solvers: List[Solver]

    # Define the simulation parameters
    df_ref_peptides: pd.DataFrame # reference database
    random_effects: bool = False
    peptide_sampling_method: Literal['', 'levenshtein', 'alanine_scanning'] = ''
    peptide_length: int = 9
    mu_immunogenic: float = 100.0
    mu_nonimmunogenic: float = 10.0
    dispersion_factor: float = 1.0
    method: Literal['threshold', 'adaptive', 'combined'] = 'threshold'
    alpha: float = 0.05 
    deconvolution_methods: List = field(default_factory=lambda: ['empirical', 'lasso', 'em'])
    min_peptide_activity: float = 10.0 # minimum to be considered a hit by 1-shot deconvolution
    _threshold: float = None

    # Define the hardware parameters
    num_processes: int = 1
        
    @property
    def threshold(self):
        return self._threshold

    @staticmethod
    def generate_random_seed() -> int:
        """
        Generates and returns a random seed.

        Returns
        -------
        random_seed :   Random seed.        
        """
        return random.randint(1, 100000000)
        
    def __post_init__(self):
        # Set threshold
        if self.random_effects:
            # If using random effects set the threshold to the
            # immunogenic mean + correction factor can get complicated
            # if num_peptides_per_pool * mu_nonimmunogenic >= mu_immunogenic,
            # so we apply a correction factor to the threshold to account for
            # the non-immunogenic peptides of n-1
            self._threshold = self.mu_immunogenic + \
                              self.mu_nonimmunogenic * max(0, self.num_peptides_per_pool-2)
        else:
            # In the non-random effects case, all the non-immunogenic peptides
            # are assigned a spot count of 0 so the threshold is set to 1
            self._threshold = 1
        assert self.num_positives <= self.num_peptides, \
            "Number of positives must be less than or equal to the number of peptides."
        assert self.num_peptides >= self.num_peptides_per_pool, \
            "Number of peptides must be greater than or equal to the number of peptides per pool."

    @staticmethod
    def evaluate_deconvolution(
            hit_peptide_ids: List[str],
            label_dict: Dict[str, int]
    ) -> Tuple[float, float, float]:
        """
        Evaluate the deconvolution results.

        Parameters
        ----------
        hit_peptide_ids     :   List of hit peptide IDs.
        label_dict          :   Dictionary of ground truth labels;
                                key     = peptide ID.
                                value   = 0 or 1.

        Returns
        -------
        sensitivity         :   Sensitivity.
        specificity         :   Specificity.
        precision           :   Precision.
        """
        positives = [peptide_id for peptide_id in label_dict.keys() if label_dict[peptide_id] == 1]
        negatives = [peptide_id for peptide_id in label_dict.keys() if label_dict[peptide_id] == 0]

        # Sensitivity: True positives / All ground truth positives aka Recall
        sensitivity = len(set(positives).intersection(set(hit_peptide_ids))) / len(positives)
        # Specificity: True negatives / All ground truth negatives aka 1 - False Positive Rate
        not_hit_peptide_ids = [peptide_id for peptide_id in label_dict.keys() if peptide_id not in hit_peptide_ids]
        specificity = len(set(not_hit_peptide_ids).intersection(set(negatives))) / len(negatives)
        # Precision: True positives / All positives predicted (PPV)
        if len(hit_peptide_ids) == 0:
            precision = 0.0
        else:
            precision = len(set(positives).intersection(set(hit_peptide_ids))) / len(hit_peptide_ids)
        return sensitivity, specificity, precision
    
    def evaluate_1shot_deconvolution(
            self,
            df_predicted_peptide_activities: pd.DataFrame,
            label_dict: Dict[str, int]
    ) -> float:
        """
        Evaluate the single-shot deconvolution results via AUCROC.

        Parameters
        ----------
        df_predicted_peptide_activities :   pd.DataFrame with the following columns:
                                            'peptide_id'
                                            'peptide_activity'
        label_dict                      :   Dictionary of ground truth labels;
                                            key     = peptide ID.
                                            value   = 0 or 1.

        Returns
        -------
        aucroc_score                   :   AUCROC score.
        """
        y_true = []
        y_predicted = []
        for peptide_id in label_dict.keys():
            peptide_activity_level = df_predicted_peptide_activities.loc[
                df_predicted_peptide_activities['peptide_id'] == peptide_id,
                'peptide_activity_level'
            ].values[0]
            y_true.append(label_dict[peptide_id])
            y_predicted.append(float(peptide_activity_level))
        return roc_auc_score(y_true, y_predicted)

    def run_worker(self, iteration: int) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Runs one iteration of the simulation. Loop through the solvers and
        generate a configuration for each solver. Then, simulate the spot
        counts for each configuration and identify the hit peptides. Finally,
        run the deconvolution on the hit peptides and return the results.

        Parameters
        ----------
        iteration       :   Number of iterations to run.

        Returns
        -------
        df_results      :   pd.DataFrame with the following columns:
                            'experiment_id'
                            'random_seed'
                            'iteration'
                            'peptides'
                            'num_peptides'
                            'num_peptides_per_pool'
                            'num_coverage'
                            'num_pools'
                            'solver'
                            'preferred_peptide_pairs'
                            'predicted_total_pools'
                            'num_violations'
                            'num_positive_peptide_sequences'
                            'positive_peptide_sequences'
                            'num_total_pools_empirical'
                            'sensitivity_empirical'
                            'specificity_empirical'
                            'precision_empirical'
                            'aucroc_score_empirical'
                            'sensitivity_lasso'
                            'specificity_lasso'
                            'precision_lasso'
                            'aucroc_score_lasso'
                            'sensitivity_em'
                            'specificity_em'
                            'precision_em'
                            'aucroc_score_em'
        df_assignments  :   pd.DataFrame with the following columns:
                            'experiment_id'
                            'iteration'
                            'solver'
                            'coverage_id'
                            'pool_id'
                            'peptide_id'
                            'peptide_sequence'
        """
        # Step 1. Generate a random seed.
        random_seed = Experiment.generate_random_seed()

        # Step 2. Sample the peptides. Return a dataframe ids and sequences and ids:labels dict
        df_peptides, label_dict = self.sample_peptides(random_seed=random_seed)
        positive_label_indices = [idx for idx, v in enumerate(label_dict.values()) if v == 1]
        df_positive_peptides = df_peptides.iloc[positive_label_indices]
        positive_peptide_sequences = df_positive_peptides['peptide_sequence']
        peptides = convert_dataframe_to_peptides(df_peptides=df_peptides)

        # Step 3. Generate a list of assignments and preferred peptides from each solver
        list_block_assignments = []
        list_preferred_peptides_pairs = []
        for solver in self.solvers:
            block_assignment, preferred_peptide_pairs = solver.generate_assignment(
                peptides=convert_dataframe_to_peptides(df_peptides=df_peptides),
                num_peptides_per_pool=self.num_peptides_per_pool,
                num_coverage=self.coverage,
                random_seed=random_seed
            )
            list_block_assignments.append(block_assignment)
            list_preferred_peptides_pairs.append(preferred_peptide_pairs)

        # Step 4. Simulate the ELISpot assay on the peptide level (configuration independent)
        peptide_spot_counts = self.simulate_peptide_spot_counts(label_dict=label_dict)

        # Step 5. Evaluate each solver and deconvolution method
        results_dict = {
            'experiment_id': [],
            'random_seed': [],
            'iteration': [],
            'peptides': [],
            'num_peptides': [],
            'num_peptides_per_pool': [],
            'num_coverage': [],
            'num_pools': [],
            'solver': [],
            'preferred_peptide_pairs': [],
            'predicted_total_pools': [],
            'num_violations': [],
            'num_positive_peptide_sequences': [],
            'positive_peptide_sequences': [],
            'num_total_pools_empirical': [],
            'sensitivity_empirical': [],
            'specificity_empirical': [],
            'precision_empirical': [],
            'aucroc_score_empirical': [],
            'sensitivity_lasso': [],
            'specificity_lasso': [],
            'precision_lasso': [],
            'aucroc_score_lasso': [],
            'sensitivity_em': [],
            'specificity_em': [],
            'precision_em': [],
            'aucroc_score_em': []
        }
        df_assignments = pd.DataFrame()
        for i in range(0, len(list_block_assignments)):
            # Index into the current block assignment and solver
            curr_block_assignment = list_block_assignments[i]
            df_assignment = list_block_assignments[i].to_dataframe()
            solver = self.solvers[i]
            preferred_peptide_pairs = list_preferred_peptides_pairs[i]
            
            # Compute the pool spot counts per pool assigned by the solver
            df_readout = Experiment.aggregate_pool_spot_counts(
                df_assignment=df_assignment,
                peptide_spot_counts=peptide_spot_counts
            )

            # Deconvolve the hit peptides
            num_total_pools_empirical = 0
            for deconvolution_method in self.deconvolution_methods:
                deconvolution_result = run_ace_deconvolve(
                    df_readout=df_readout,
                    block_assignment=curr_block_assignment,
                    mode=deconvolution_method,
                    statistical_min_peptide_activity=self.min_peptide_activity,
                    empirical_min_coverage=1,
                    empirical_min_spot_count=int(self._threshold),
                    verbose=False
                )
                df_deconvolution_result = deconvolution_result.to_dataframe()

                # Evaluate the deconvolution results
                aucroc_score = self.evaluate_1shot_deconvolution(
                    df_predicted_peptide_activities=df_deconvolution_result,
                    label_dict=label_dict
                )

                # Filter empirical deconvolution for self.coverage peptide activity level
                if deconvolution_method == DeconvolveModes.EMPIRICAL:
                    df_deconvolution_result = df_deconvolution_result.loc[
                        df_deconvolution_result['peptide_activity_level'] == self.coverage,:
                    ]
                    candidate_peptide_ids = df_deconvolution_result.loc[
                        df_deconvolution_result['deconvolution_result'] == DeconvolutionLabels.CANDIDATE_HIT,
                        'peptide_id'
                    ].values.tolist()
                    num_total_pools_empirical = len(df_assignment['pool_id'].unique()) + len(candidate_peptide_ids)

                hit_peptide_ids = df_deconvolution_result.loc[
                    df_deconvolution_result['deconvolution_result'].isin(
                        [DeconvolutionLabels.CONFIDENT_HIT, DeconvolutionLabels.CANDIDATE_HIT]
                    ),
                    'peptide_id'
                ].values.tolist()
                sensitivity, specificity, precision = Experiment.evaluate_deconvolution(
                    hit_peptide_ids=hit_peptide_ids,
                    label_dict=label_dict
                )
                results_dict['sensitivity_%s' % deconvolution_method].append(sensitivity)
                results_dict['specificity_%s' % deconvolution_method].append(specificity)
                results_dict['precision_%s' % deconvolution_method].append(precision)
                results_dict['aucroc_score_%s' % deconvolution_method].append(aucroc_score)

            # Append results to results_dict
            results_dict['experiment_id'].append(self.experiment_id)
            results_dict['random_seed'].append(random_seed)
            results_dict['iteration'].append(iteration)
            results_dict['solver'].append(solver.name)
            results_dict['peptides'].append(';'.join(['%s,%s' % (peptide_id, peptide_sequence) for peptide_id, peptide_sequence in peptides]))
            results_dict['num_peptides'].append(len(peptides))
            results_dict['num_peptides_per_pool'].append(self.num_peptides_per_pool)
            results_dict['num_coverage'].append(self.coverage)
            results_dict['num_pools'].append(len(df_assignment['pool_id'].unique()))
            results_dict['preferred_peptide_pairs'].append(';'.join(['%s,%s' % (peptide_id, peptide_sequence) for peptide_id, peptide_sequence in preferred_peptide_pairs]))
            results_dict['predicted_total_pools'].append(num_total_pools_empirical)
            results_dict['num_violations'].append(curr_block_assignment.num_violations())
            results_dict['num_positive_peptide_sequences'].append(len(positive_peptide_sequences))
            results_dict['positive_peptide_sequences'].append(';'.join(positive_peptide_sequences))
            results_dict['num_total_pools_empirical'].append(num_total_pools_empirical)

            # Concatenate assignments
            df_assignment['experiment_id'] = [self.experiment_id] * len(df_assignment)
            df_assignment['iteration'] = [iteration] * len(df_assignment)
            df_assignment['solver'] = [solver.name] * len(df_assignment)
            df_assignments = pd.concat([df_assignments, df_assignment])
        df_results = pd.DataFrame(results_dict)
        return df_results, df_assignments
        
    def run(self, num_iterations: int = 1) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Runs an ELISpot simulation.

        Parameters
        ----------
        num_iterations      :   Number of iterations (default: 1).

        Returns
        -------
        df_results_all      :   pd.DataFrame with the following columns:
                                'experiment_id'
                                'random_seed'
                                'iteration'
                                'peptides'
                                'num_peptides'
                                'num_peptides_per_pool'
                                'num_coverage'
                                'num_pools'
                                'solver'
                                'preferred_peptide_pairs'
                                'predicted_total_pools'
                                'num_violations'
                                'num_positive_peptide_sequences'
                                'positive_peptide_sequences'
                                'num_total_pools_empirical'
                                'sensitivity_empirical'
                                'specificity_empirical'
                                'precision_empirical'
                                'aucroc_score_empirical'
                                'sensitivity_lasso'
                                'specificity_lasso'
                                'precision_lasso'
                                'aucroc_score_lasso'
                                'sensitivity_em'
                                'specificity_em'
                                'precision_em'
                                'aucroc_score_em'
        df_assignments_all  :   pd.DataFrame with the following columns:
                                'experiment_id'
                                'iteration'
                                'solver'
                                'coverage_id'
                                'pool_id'
                                'peptide_id'
                                'peptide_sequence'
        """
        iterations = [i for i in range(1, num_iterations + 1)]

        # Multiprocess experiments k iterations
        pool = mp.Pool(processes=self.num_processes)
        func = partial(self.run_worker)
        results = pool.map(func, iterations)
        pool.close()
        
        # Postprocess results
        df_results_all = pd.DataFrame()
        df_assignments_all = pd.DataFrame()
        for df_results, df_assignments in results:
            df_results_all = pd.concat([df_results_all, df_results])
            df_assignments_all = pd.concat([df_assignments_all, df_assignments])
        return df_results_all, df_assignments_all

    def sample_peptides(self, random_seed: int) -> Tuple[pd.DataFrame, Dict[str, int]]:
        """
        Sample the peptide sequences and their immunogenicity status
        from the reference data.

        NOTE: We do this instead of randomly initializing peptide sequences
        because the ACE Neural Engine was trained on real sequence data. Thus,
        to make meaningful disallowed peptide pairings we utilize a reference
        dataset.

        Parameters
        ----------
        random_seed :   Random seed.

        Returns
        -------
        df_peptides :   pd.DataFrame with the following columns:
                        'peptide_id'
                        'peptide_sequence'
        labels_dict :   Dictionary;
                        key     = peptide ID.
                        value   = 0 or 1.
        """
        if self.peptide_sampling_method == '':
            # Sample the peptides from the reference dataframe without replacement
            pos_df = self.df_ref_peptides[self.df_ref_peptides['Binding'] == 1].sample(n=self.num_positives, random_state=random_seed)
            neg_df = self.df_ref_peptides[self.df_ref_peptides['Binding'] == 0].sample(n=self.num_peptides - self.num_positives, random_state=random_seed)
        elif self.peptide_sampling_method == 'levenshtein':
            peptides_df = copy.deepcopy(self.df_ref_peptides)
            # Pop a random row from the reference dataframe
            sampled_peptide = peptides_df.sample(n=1, random_state=random_seed)
            peptides_df = peptides_df.drop(sampled_peptide.index)
            # Sort peptides on distance to the sampled peptide
            peptides_df['distance'] = peptides_df['Epitope'].apply(lambda x: Experiment.levenshtein_distance(x, sampled_peptide['Epitope'].values[0]))
            peptides_df = peptides_df.sort_values(by='distance', inplace=False)
            # Select the top n peptides from the sorted dataframe stratified by immunogenicity
            pos_df = peptides_df[peptides_df['Binding'] == 1].head(self.num_positives)
            neg_df = peptides_df[peptides_df['Binding'] == 0].head(self.num_peptides - self.num_positives)
        elif self.peptide_sampling_method == 'alanine_scanning':
            length_filter_mask = (self.df_ref_peptides['Epitope'].str.len() == self.peptide_length)
            length_filtered_df_ref_peptides = self.df_ref_peptides.loc[length_filter_mask]
            num_pos_ref_peptides = math.ceil(self.num_positives / self.peptide_length)
            num_neg_ref_peptides = int(self.num_peptides / self.peptide_length) - num_pos_ref_peptides
            pos_df = length_filtered_df_ref_peptides[length_filtered_df_ref_peptides['Binding'] == 1].sample(n=num_pos_ref_peptides, random_state=random_seed)
            neg_df = length_filtered_df_ref_peptides[length_filtered_df_ref_peptides['Binding'] == 0].sample(n=num_neg_ref_peptides, random_state=random_seed)
            pos_data = {
                'Epitope': [],
                'Binding': []
            }
            num_positives_ = self.num_positives
            for index, row in pos_df.iterrows():
                peptide_sequence = row['Epitope']
                peptide_sequences = Experiment.perform_alanine_scanning(peptide_sequence=peptide_sequence)
                if num_positives_ > 0:
                    if num_positives_ >= len(peptide_sequences):
                        for peptide_sequence_ in peptide_sequences:
                            pos_data['Epitope'].append(peptide_sequence_)
                            pos_data['Binding'].append(1)
                            num_positives_ -= 1
                    else :
                        positive_indices = random.sample(list(range(0, len(peptide_sequences))), num_positives_)
                        for idx in range(0, len(peptide_sequences)):
                            peptide_sequence_ = peptide_sequences[idx]
                            pos_data['Epitope'].append(peptide_sequence_)
                            if idx in positive_indices:
                                pos_data['Binding'].append(1)
                                num_positives_ -= 1
                            else:
                                pos_data['Binding'].append(0)
                else:
                    for peptide_sequence_ in peptide_sequences:
                        pos_data['Epitope'].append(peptide_sequence_)
                        pos_data['Binding'].append(0)
            pos_df = pd.DataFrame(pos_data)
            neg_data = {
                'Epitope': [],
                'Binding': []
            }
            for index, row in neg_df.iterrows():
                peptide_sequence = row['Epitope']
                peptide_sequences = Experiment.perform_alanine_scanning(peptide_sequence=peptide_sequence)
                for peptide_sequence_ in peptide_sequences:
                    neg_data['Epitope'].append(peptide_sequence_)
                    neg_data['Binding'].append(0)
            neg_df = pd.DataFrame(neg_data)

        # Concatenate the positive and negative peptides and shuffle the dataframe
        peptides_df = pd.concat([pos_df, neg_df], axis=0).sample(frac=1.0, random_state=random_seed)
        # Reset the index of the dataframe
        peptides_df = peptides_df.reset_index(drop=True)
        # Add a peptide_id column to the dataframe
        peptides_df['peptide_id'] = [f'peptide_{i}' for i in range(1, len(peptides_df) + 1)]
        
        # Select the peptide sequences and peptide ids
        df_peptides = peptides_df[['peptide_id', 'Epitope']]
        df_peptides.columns = ['peptide_id', 'peptide_sequence']

        # Create a dictionary of labels
        labels_dict = {peptide_id: label for peptide_id, label in zip(peptides_df['peptide_id'], peptides_df['Binding'])}
        return df_peptides, labels_dict
    
    @staticmethod
    def levenshtein_distance(sequence_1: str, sequence_2: str) -> int:
        """
        Calculates the Levenshtien distance between two sequences.

        Parameters
        ----------
        sequence_1      :   Sequence.
        sequence_2      :   Sequence.

        Returns
        -------
        levenshtein     :   Levenshtein distance.
        """
        return Levenshtein.distance(sequence_1, sequence_2)
    
    @staticmethod
    def sample_spot_counts(
        mean: int, 
        dispersion_factor: float, 
        num_samples: int
    ) -> List[int]:
        """
        Sample the number of spots for a peptide pool using a negative binomial distribution.

        Derivation of NegBinom Parameters:

            Let X be the number RV of spots sampled for a peptide pool.
            Let p be the probability of sampling a spot for a peptide in the pool.
            Let r be the number of successes (spots) we want to sample.
            Let k be the number of failures (non-spots) we want to sample.

            Then, the probability mass function for X is given by:

            P(X=k) = (k+r-1)C(k) * p^r * (1-p)^k
            
            where (k+r-1)C(k) is the binomial coefficient.

            The mean and variance of X are given by:

            mean = r(1-p)/p
            var = r(1-p)/p^2

            We can solve for p and r in terms of the mean and variance:

            p = mean / var
            r = mean^2 / (var - mean)
        
        Parameters
        ----------
        mean                :   Mean sample spot count.
        dispersion_factor   :   Dispersion factor.
        num_samples         :   Number of samples.

        Returns
        -------
        spot_counts         :   List of integers. The number of elements in this
                                list is the same as 'num_samples'.
        """
        if dispersion_factor == 0:        
            return [mean]*num_samples
        elif dispersion_factor < 1:
            raise ValueError("Dispersion_factor must be an integer greater than or equal to 1, or 0.")
        elif dispersion_factor == 1:
            return list(np.random.poisson(mean, num_samples))
        else:    
            variance = mean*dispersion_factor
            p = mean/variance
            r = mean**2/(variance-mean)
            return list(np.random.negative_binomial(r, p, num_samples))
    
    @staticmethod
    def perform_alanine_scanning(peptide_sequence: str) -> List[str]:
        """
        Performs alanine scanning and returns a list of peptide sequences.

        Parameters
        ----------
        peptide_sequence    :   Peptide sequence.

        Returns
        -------
        peptide_sequences   :   List of peptide sequences.
        """
        peptide_sequences = []
        for i in range(0, len(peptide_sequence)):
            peptide_sequence_ = list(peptide_sequence)
            peptide_sequence_[i] = 'A'
            peptide_sequences.append(''.join(peptide_sequence_))
        return peptide_sequences

    def simulate_peptide_spot_counts(
            self,
            label_dict: Dict[str, int]
    ) -> Dict[str, List[int]]:
        """
        Given a list of peptide_ids and their immunogenicity labels,
        simulate the number of spots sampled for each peptide for all
        coverages.

        Parameters
        ----------
        label_dict      :   Dictionary of ground truth labels;
                            key     = peptide ID.
                            value   = label (0 or 1).

        Returns
        -------
        spot_counts     :   Dictionary;
                            key     = peptide ID.
                            value   = List[int];
                                      each element corresponds to the spot count
                                      at the indexed coverage.
        """
        # Step 1. Create dictionary to store peptide spot count per coverage
        # key   = peptide ID
        # value = list of spot counts (per coverage)
        peptide_spot_counts = {}

        # Step 2. Determine immunogenic and non-immunogenic mean and dispersion factor values
        if self.random_effects:
            # Simulate using random_effects (negative binom)
            immunogenic_mean = self.mu_immunogenic
            non_immunogenic_mean = self.mu_nonimmunogenic
            dispersion_factor = self.dispersion_factor
        else:
            # Simulate using deterministic (binary) parameters
            immunogenic_mean = 1
            non_immunogenic_mean = 0
            dispersion_factor = 0

        # Step 3. Simulate peptide spot counts
        for peptide_id in label_dict.keys():
            if label_dict[peptide_id] == 1:
                peptide_spot_counts[peptide_id] = Experiment.sample_spot_counts(
                    mean=immunogenic_mean,
                    dispersion_factor=dispersion_factor,
                    num_samples=self.coverage
                )
            else:
                peptide_spot_counts[peptide_id] = Experiment.sample_spot_counts(
                    mean=non_immunogenic_mean,
                    dispersion_factor=dispersion_factor,
                    num_samples=self.coverage
                )
        return peptide_spot_counts
    
    @staticmethod
    def aggregate_pool_spot_counts(
            df_assignment: pd.DataFrame,
            peptide_spot_counts: Dict[str, List[int]]
    ) -> pd.DataFrame:
        """
        Aggregate the ELISPOT assay using the optimal peptide pools configuration.
        Determine how spots are sampled for each peptide pool. 

        Assumptions:
        1.  The number of spots sampled for each peptide pool follows a
            negative binomial distribution.
        2.  The number of spots sampled for each peptide pool is independent
            of the other peptide pools.
        3.  Immunogenic and non-immunogenic peptides have the same shape
            parameters for the negative binomial distribution, but different
            mean parameters. This is because we assume the experimental
            error/noise is the same for both whereas true immunogenic
            peptides will likely have more spots than non-immunogenic peptides.
        4.  We assume that having multiple immunogenic peptides in a pool
            additively increases the number of spots sampled.

        Parameters
        ----------
        df_assignment       :   pd.DataFrame with the following columns:
                                'coverage_id'
                                'pool_id'
                                'peptide_id'
        peptide_spot_counts :   Dictionary;
                                key     =   peptide ID.
                                value   =   List[int]; each element corresponds
                                            to the spot count at the indexed
                                            coverage.

        Returns
        -------
        df_readout          :   pd.DataFrame with the following columns:
                                'pool_id'
                                'spot_count'
        """
        # Step 1. Aggregate peptide spot count to pool spot counts
        pool_spot_counts = {} # key = pool ID, value = spot count
        for peptide_id in peptide_spot_counts.keys():
            pool_ids = df_assignment.loc[
                df_assignment['peptide_id'] == peptide_id,
                'pool_id'
            ].values.tolist()
            curr_peptide_spot_counts = peptide_spot_counts[peptide_id]
            random.shuffle(pool_ids)
            for i in range(0, len(pool_ids)):
                curr_pool_id = pool_ids[i]
                curr_peptide_spot_count = curr_peptide_spot_counts[i]
                if curr_pool_id not in pool_spot_counts.keys():
                    pool_spot_counts[curr_pool_id] = 0
                pool_spot_counts[curr_pool_id] += curr_peptide_spot_count

        # Step 2. Dump data into DataFrame
        data = {
            'pool_id': [],
            'spot_count': []
        }
        for key, val in pool_spot_counts.items():
            data['pool_id'].append(key)
            data['spot_count'].append(val)
        return pd.DataFrame(data)

    def deconvolve_hit_pools(
            self, 
            pool_spot_counts: Dict[int, int], 
            method: Literal['threshold', 'adaptive', 'combined'],
            threshold: float,
            alpha: float = 0.05,
            negative_control_sampling_iters: int = 10000
    ) -> List[int]:
        """
        Identify the hit pools given the simulated spot counts for each pool.
        Pass the hit pools to the deconvolution method implemented by the
        different solvers. Based on the Empirical Rule (ER) response definition
        criteria from hte Moodie et al. (2010) paper: Response definition
        criteria for ELISPOT assays revisited:
        https://pubmed.ncbi.nlm.nih.gov/20549207/

        Parameters
        ----------
        pool_spot_counts    :   Dictionary;
                                key     = pool ID.
                                value   = total spot count.
        method              :   Method (allowed values: 'threshold', 'adaptive', 'combined').

                                If 'threshold', hit pools are those with
                                spot counts greater than or equal to the threshold
                                (Dubey et al.).

                                If 'adaptive', hit pools are identified based
                                on the control pools (e.g. DMSO).

                                If 'combined', hit pools are identified based
                                on the control pools (e.g. DMSO) and the threshold.
        threshold           :   Threshold. If unspecified, the threshold is set
                                by default based on whether random effects are
                                taken into account.
        alpha               :   Significance level (default: 0.05).

        Returns
        -------
        hit_pool_ids        :   List of hit pool IDs.
        """
        significance = (1 - alpha) * 100
        if method == 'threshold':
            hit_pool_ids = []
            for pool_id, spot_count in pool_spot_counts.items():
                if spot_count >= threshold:
                    hit_pool_ids.append(pool_id)
        elif method == 'adaptive':
            assert 0 < alpha < 1, "alpha must be between 0 and 1."
            if not self.random_effects:
                raise ValueError("Adaptive thresholding only works with random effects")
                exit(1)
            # Simulate modified DMSO controls
            # Normal DMSO controls are simulated with the non-immunogenic mean across the whole pool
            # Here we simulate the DMSO controls with the non-immunogenic mean per peptide and add across the whole pool
            negative_controls = []
            for _ in range(0, negative_control_sampling_iters):
                negative_controls.append(sum(
                    Experiment.sample_spot_counts(
                        mean=self.mu_nonimmunogenic,
                        dispersion_factor=self.dispersion_factor,
                        num_samples=self.num_peptides_per_pool
                )))
            hit_pool_ids = []
            for pool_id, spot_count in pool_spot_counts.items():
                if spot_count >= np.percentile(negative_controls, significance):
                    hit_pool_ids.append(pool_id)
        elif method == 'combined':
            thresholded_pool_ids = []
            for pool_id, spot_count in pool_spot_counts.items():
                if spot_count >= threshold:
                    thresholded_pool_ids.append(pool_id)
            negative_controls = []
            for _ in range(negative_control_sampling_iters):
                negative_controls.append(sum(
                    Experiment.sample_spot_counts(
                        mean=self.mu_nonimmunogenic,
                        dispersion_factor=self.dispersion_factor,
                        num_samples=self.num_peptides_per_pool
                )))
            adaptive_pool_ids = []
            for pool_id, spot_count in pool_spot_counts.items():
                if spot_count >= np.percentile(negative_controls, significance):
                    adaptive_pool_ids.append(pool_id)
            hit_pool_ids = list(set(thresholded_pool_ids + adaptive_pool_ids))
        else:
            raise ValueError("Invalid method. Please choose from ['threshold', 'adaptive', or 'combined'].")
            exit(1)
        return hit_pool_ids

