"""
The purpose of this python3 script is to implement the Experiment dataclas
"""


import copy
import multiprocessing as mp
import numpy as np
import pandas as pd
import random
import time
from dataclasses import dataclass, field
from functools import partial
from typing import Dict, List, Tuple
from acelib.main import run_ace_deconvolve
from .solver import Solver
import Levenshtein


@dataclass
class Experiment:
    """
    Class to simulate an end-to-end ELISpot experiment. Class accepts
    a Solver objects, which are used to generate a configuration
    for a given sample of peptides. For each configuration, the class
    simulates the number of spots per well in the experiment. Then, the 
    putative positive peptides are identified and test statistics are
    calculated for the deconvolution results.

    Simulation Procedure:
    1. Initialize the object with the given assay parameters.
    2. Load the reference database of peptides to sample from
    3. Generate a configuration of peptides to be used in the experiment
    4. Simulate the spot counts for each well in the experiment
    5. Identify the hit peptides from the simulated spot counts
    6. Use the candidate + hit peptides list to run the deconvolution
    7. Return statistics on the deconvolution results
    """
    
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
    peptide_scan: bool = False
    mu_immunogenic: float = 100.0
    mu_nonimmunogenic: float = 10.0
    dispersion_factor: float = 1.0
    method: str = 'threshold' # Allowed values: 'threshold', 'adaptive', 'combined'
    alpha: float = 0.05
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
        return random.randint(1, 1000000)
        
    def __post_init__(self):
        # Set threshold
        if self.random_effects:
            # If using random effects set the threshold to the immunogenic mean
            # Can get complicated if num_peptides_per_pool * mu_nonimmunogenic >= mu_immunogenic so we 
            # apply a correction factor to the threshold to account for the non-immunogenic peptides of n-1
            self._threshold = self.mu_immunogenic + self.mu_nonimmunogenic*max(0, self.num_peptides_per_pool-2)
        else:
            # In the non-random effects case, all the non-immunogenic peptides 
            # are assigned a spot count of 0 so the threshold is set to 1
            self._threshold = 1

        assert self.num_positives <= self.num_peptides, "Number of positives must be less than or equal to the number of peptides."
        assert self.num_peptides >= self.num_peptides_per_pool, "Number of peptides must be greater than or equal to the number of peptides per pool."

    def evaluate_deconvolution(
            self, 
            hit_peptide_ids: List[str],
            label_dict: Dict[str, int]
    ) -> Tuple[float, float, float]:
        """
        Evaluate the deconvolution results.

        Parameters
        ----------
        hit_peptide_ids         :   List of hit peptide IDs.
        label_dict      :   Dictionary of ground truth labels;
        """
        positives = [peptide_id for peptide_id in label_dict.keys() if label_dict[peptide_id] == 1]
        negatives = [peptide_id for peptide_id in label_dict.keys() if label_dict[peptide_id] == 0]
        # Sensitivity: True positives / All ground truth positives aka Recall
        sensitivity = len(set(positives).intersection(set(hit_peptide_ids))) / len(positives)
        # Specificity: True negatives / All ground truth negatives aka 1 - False Positive Rate
        not_hit_peptide_ids = [peptide_id for peptide_id in label_dict.keys() if peptide_id not in hit_peptide_ids]
        specificity = len(set(not_hit_peptide_ids).intersection(set(negatives))) / len(negatives)
        # Precision: True positives / All positives predicted (PPV)
        precision = len(set(positives).intersection(set(hit_peptide_ids))) / len(hit_peptide_ids)
        return sensitivity, specificity, precision

    def run_worker(self, iteration: int) -> dict:
        """
        Run one iteration of the simulation. Loop through the solvers and
        generate a configuration for each solver. Then, simulate the spot
        counts for each configuration and identify the hit peptides. Finally,
        run the deconvolution on the hit peptides and return the results.
        """
        # Step 1. Generate a random seed.
        random_seed = Experiment.generate_random_seed()
        
        # Step 2. Sample the peptides. Return a dataframe ids and sequences and ids:labels dict
        df_peptides, label_dict = self.sample_peptides(random_seed=random_seed)

        # Step 3. Generate a configuration for each solver
        list_configurations = []
        for solver in self.solvers:
            print("Starting SOLVER: %s" % solver.name)
            df_configuration = solver.generate_configuration(
                df_peptides=df_peptides,
                num_peptides_per_pool=self.num_peptides_per_pool,
                num_coverage=self.coverage
            )
            list_configurations.append(df_configuration)
            print("Finished SOLVER: %s" % solver.name)

        # Step 4. Simulate the ELISPOT assay on the Peptide Level (Configuration Independent)
        peptide_spot_counts = self.simulate_peptide_spot_counts(label_dict=label_dict)

        # Step 5. Evaluate each solver
        results_data = {
            'iteration': [],
            'solver': [],
            'num_total_pools': [],
            'sensitivity': [],
            'specificity': [],
            'precision': []
        }
        for i in range(0, len(list_configurations)):
            df_configuration = list_configurations[i]
            solver = self.solvers[i]
            pool_spot_counts = self.aggregate_pool_spot_counts(
                df_configuration=df_configuration, 
                peptide_spot_counts=peptide_spot_counts
            )
            hit_pool_ids = self.deconvolve_hit_pools(
                pool_spot_counts=pool_spot_counts,
                method=self.method,
                threshold=self.threshold,
                alpha=self.alpha
            )
            df_hits = self.deconvolve_hit_peptides(
                hit_pool_ids=hit_pool_ids,
                df_configuration=df_configuration,
                min_coverage=self.coverage
            )
            candidate_peptide_ids = df_hits.loc[df_hits['deconvolution_result'] == 'candidate_hit','peptide_id'].values.tolist()
            num_total_pools_solver = len(df_configuration['pool_id'].unique()) + len(candidate_peptide_ids)
            sensitivity, specificity, precision = self.evaluate_deconvolution(
                hit_peptide_ids=df_hits['peptide_id'].values.tolist(),
                label_dict=label_dict
            )

            # Append results to results_data dict for that solver
            results_data['iteration'].append(iteration)
            results_data['solver'].append(solver.name)
            results_data['num_total_pools'].append(num_total_pools_solver)
            results_data['sensitivity'].append(sensitivity)
            results_data['specificity'].append(specificity)
            results_data['precision'].append(precision)

        return results_data
        
    def run(self, num_iterations=1) -> pd.DataFrame:
        """
        Runs an ELISpot simulation.

        Parameters
        ----------
        num_iterations  :   Number of iterations (default: 1).

        Returns
        -------
        df_results      :   pd.DataFrame with the following columns:
                            '
        """
        iterations = [i for i in range(1, num_iterations + 1)]

        # Multiprocess experiments k iterations
        pool = mp.Pool(processes=self.num_processes)
        func = partial(self.run_worker)
        results = pool.map(func, iterations)
        pool.close()
        
        # Postprocess results
        df_results = pd.DataFrame(results)
        return df_results

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
                        key = peptide ID.
                        value = label (0 or 1).
        """
        if self.peptide_scan:
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
        else:
            # Sample the peptides from the reference dataframe
            pos_df = self.df_ref_peptides[self.df_ref_peptides['Binding'] == 1].sample(n=self.num_positives, random_state=random_seed)  
            neg_df = self.df_ref_peptides[self.df_ref_peptides['Binding'] == 0].sample(n=self.num_peptides - self.num_positives, random_state=random_seed)
        # Concatenate the positive and negative peptides and shuffle the dataframe
        peptides_df = pd.concat([pos_df, neg_df], axis=0).sample(frac=1.0, random_state=random_seed)
        # Reset the index of the dataframe
        peptides_df = peptides_df.reset_index(drop=True)
        # Add a peptide_id column to the dataframe
        peptides_df['peptide_id'] = [f'peptide_{i}' for i in peptides_df.index]
        
        # Select the peptide sequences and peptide ids
        df_peptides = peptides_df[['peptide_id', 'Epitope']]
        df_peptides.columns = ['peptide_id', 'peptide_sequence']

        # Create a dictionary of labels
        labels_dict = {peptide_id: label for peptide_id, label in zip(peptides_df['peptide_id'], peptides_df['Binding'])}
        return df_peptides, labels_dict
    
    @staticmethod
    def levenshtein_distance(seq1, seq2) -> int:
        """
        Calculates the Levenshtien distance between two sequences.

        Parameters
        ----------
        seq1            :   Sequence.
        seq2            :   Sequence.

        Returns
        -------
        levenshtein     :   Levenshtein distance.
        """
        return Levenshtein.distance(seq1, seq2)
    
    def generate_configuration(
            self, 
            solver: Solver, 
            df_peptides: pd.DataFrame,
            **kwargs
    ) -> pd.DataFrame:
        """
        Calls the solver's generate_configuration method to generate a configuration
        of peptide pools to be used in the experiment.

        Parameters
        ----------
        solver              :   Solver object.
        df_peptides         :   pd.DataFrame with the following columns:
                                'peptide_id'
                                'peptide_sequence'

        Returns
        -------
        df_configuration    :   pd.DataFrame with the following columns:
                                'coverage_id'
                                'pool_id'
                                'peptide_id'
        """
        return solver.generate_configuration(
            df_peptides=df_peptides,
            num_peptides_per_pool=self.num_peptides_per_pool,
            num_coverage=self.coverage,
            kwargs=kwargs
        )
    
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
        spot_counts         :   List of integers; 'num_samples' number of spot counts.
        """
        if dispersion_factor == 0:        
            return [mean]*num_samples
        elif dispersion_factor < 1:
            raise ValueError("dispersion_factor must be greater than or equal to 1")
        elif dispersion_factor == 1:
            return list(np.random.poisson(mean, num_samples))
        else:    
            variance = mean*dispersion_factor
            p = mean/variance
            r = mean**2/(variance-mean)
            return list(np.random.negative_binomial(r, p, num_samples))
    
    def simulate_peptide_spot_counts(
            self, 
            label_dict: dict
    ) -> dict:
        """
        Given a list of peptide sequences and their immunogenicity labels,
        simulate the number of spots sampled for each peptide for a given
        coverage.

        Parameters
        ----------
        labels          :   Dictionary of ground truth labels;
                            key = peptide ID.
                            value = label (0 or 1).

        Returns
        -------
        spot_counts     :   Dictionary; 
                            key = peptide ID.
                            value = List[int]; each element corresponds to the spot count at the indexed coverage.
        """
        peptide_ids = label_dict.keys()
        # Peptide spot counts per coverage
        peptide_spot_counts = {peptide_id:[] for peptide_id in peptide_ids}

        # Simulate using random_effects (negative binom)
        if self.random_effects:
            imm_mean = self.mu_immunogenic
            nonimm_mean = self.mu_nonimmunogenic
            disp = self.dispersion_factor
        # Simulate using deterministic (binary) parameters
        else:
            imm_mean = 1
            nonimm_mean = 0
            disp = 0
            
        for peptide_id in peptide_ids:
            if label_dict[peptide_id] == 1:
                peptide_spot_counts[peptide_id].append(Experiment.sample_spot_counts(imm_mean, disp, num_samples=self.coverage))   
            else:
                peptide_spot_counts[peptide_id].append(Experiment.sample_spot_counts(nonimm_mean, disp, num_samples=self.coverage))
        return peptide_spot_counts
    
    def aggregate_pool_spot_counts(
            self, 
            df_configuration: pd.DataFrame,
            peptide_spot_counts: Dict[str,List[int]]
    ) -> Dict[str,int]:
        """
        Aggregate the ELISPOT assay using the optimal peptide pools configuration.
        Determine how spots are sampled for each peptide pool. 

        Assumptions:
            1. The number of spots sampled for each peptide pool follows a negative binomial distribution.
            2. The number of spots sampled for each peptide pool is independent of the other peptide pools.
            3. Immunogenic and non-immunogenic peptides have the same shape parameters for the negative binomial distribution, 
            but different mean parameters. This is because we assume the experimental error/noise is the same for both whereas
            true immunogenic peptides will likely have more spots than non-immunogenic peptides.
            4. We assume that having multiple immunogenic peptides in a pool addiditively increases the number of spots sampled

        Parameters
        ----------
        df_configuration   :   pd.DataFrame with the following columns:
                                'coverage_id'
                                'pool_id'
                                'peptide_id'
        peptide_spot_counts :   Dictionary;
                                key = peptide ID.
                                value = List[int]; each element corresponds to the spot count at the indexed coverage.
                                
        Returns
        -------
        pool_spot_counts    :   Dictionary;
                                key = pool ID.
                                value = total spot count.
        """
        # Step 1. Prepare pool ID dictionary
        #pools = {}
        # TODO: Fix the peptide spot count to get the different coverage
        #for pool_id in df_configuration['pool_id'].unique():                    

        peptide_spot_counts = copy.deepcopy(peptide_spot_counts)
        pool_spot_counts = {i:0 for i in df_configuration['pool_id'].unique()}
        
        for pool_id in df_configuration['pool_id'].unique():
            pool_df = df_configuration[df_configuration['pool_id']==pool_id]
            for peptide_id in pool_df['peptide_id']:
                pool_spot_counts[pool_id] += peptide_spot_counts[peptide_id][0].pop(-1)
                
        return pool_spot_counts
    
    def deconvolve_hit_pools(
            self, 
            pool_spot_counts: Dict[str, int], 
            method,
            threshold,
            alpha=0.05
    ) -> List[str]:
        """
        Identify the hit pools given the simulated spot counts for each pool. Pass the hit pools
        to the deconvolution method implemented by the different solvers. Based on the Empirical Rule (ER) response 
        definition criteria from hte Moodie et al. (2010) paper: Response definition criteria for ELISPOT assays revisited
        https://pubmed.ncbi.nlm.nih.gov/20549207/

        Parameters
        ----------
        pool_spot_counts    :   Dictionary;
                                key = pool ID.
                                value = total spot count.
        method              :   Method (allowed values: 'threshold', 'adaptive').
                                'threshold': Identify the hit pools as those with a spot count greater than or equal to the threshold (Dubey et al.).
                                'adaptive': Identify the hit pools based on the mock wells (DMSO) controls.
                                'combined': Identify the hit pools based on the mock wells (DMSO) controls and the threshold.                
        threshold           :   Threshold (default: None). If None, the threshold is set by default based on whether or not random effects are taken into account.
        alpha               :   Significance level (default: 0.05).  
        """
        if method == 'threshold':
            hit_pools = [pool_id for pool_id, spot_count in pool_spot_counts.items() if spot_count >= threshold]
        elif method == 'adaptive':
            assert alpha > 0 and alpha < 1, "alpha must be between 0 and 1."
            significance = (1 - alpha)*100
            if not self.random_effects:
                raise ValueError("Adaptive thresholding only works with random effects")
            # Simulate modified DMSO controls
            # Normal DMSO controls are simulated with the non-immunogenic mean across the whole pool
            # Here we simulate the DMSO controls with the non-immunogenic mean per peptide and add across the whole pool
            negative_controls = [sum(Experiment.sample_spot_counts(mean=self.mu_nonimmunogenic, dispersion_factor=self.dispersion_factor, num_samples=self.num_peptides_per_pool)) for _ in range(10000)]
            hit_pools = [pool_id for pool_id, spot_count in pool_spot_counts.items() if spot_count >= np.percentile(negative_controls, significance)]
        elif method == 'combined':
            thresholded_pools = [pool_id for pool_id, spot_count in pool_spot_counts.items() if spot_count >= threshold]
            negative_controls = [sum(Experiment.sample_spot_counts(mean=self.mu_nonimmunogenic, dispersion_factor=self.dispersion_factor, num_samples=self.num_peptides_per_pool)) for _ in range(10000)]
            adaptive_pools = [pool_id for pool_id, spot_count in pool_spot_counts.items() if spot_count >= np.percentile(negative_controls, significance)]
            hit_pools = list(set(thresholded_pools + adaptive_pools))
        else:
            raise ValueError("Invalid method. Please choose from ['threshold', 'adaptive', or 'combined'].")
        return hit_pools

    def deconvolve_hit_peptides(
            self,
            hit_pool_ids: List[str],
            df_configuration: pd.DataFrame,
            min_coverage: int
    ) -> pd.DataFrame:
        """
        Deconvolves hit peptide IDs.

        Parameters
        ----------
        hit_pool_ids        :   List of pool IDs.
        df_configuration    :   pd.DataFrame with the following columns:
                                'coverage_id'
                                'pool_id'
                                'peptide_id'
                                'peptide_sequence'
        min_coverage        :   Minimum coverage.
    
        Returns
        -------
        df_hits             :   DataFrame with the following columns:
                                'peptide_id'
                                'peptide_sequence'
                                'pool_ids'
                                'num_coverage'
                                'deconvolution_result'
        """
        return run_ace_deconvolve(
            hit_pool_ids=hit_pool_ids,
            df_configuration=df_configuration,
            min_coverage=min_coverage
        )
        
    