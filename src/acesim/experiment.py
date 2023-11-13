"""
The purpose of this python3 script is to implement the Experiment dataclass.
"""


import Levenshtein
import math
import numpy as np
import pandas as pd
import random
from typing import List, Literal, Tuple
from dataclasses import dataclass
from sklearn.metrics import roc_auc_score
from acelib.block_assignment import BlockAssignment
from acelib.deconvolution import DeconvolutionResult
from acelib.constants import DeconvolutionLabels
from acesim.utilities import *


@dataclass
class Experiment:

    @staticmethod
    def generate_random_seed() -> int:
        """
        Generates and returns a random seed.

        Returns
        -------
        random_seed :   Random seed.
        """
        return random.randint(1, 100000000)

    @staticmethod
    def sample_peptides(
            df_ref_peptides: pd.DataFrame,
            num_peptides: int,
            num_peptides_immunogenic: int,
            peptide_length: int,
            sampling_method: Literal['random', 'levenshtein', 'alanine_scanning']
    ) -> pd.DataFrame:
        """
        Samples peptides.

        Parameters
        ----------
        df_ref_peptides             :   Pandas DataFrame of reference peptides
                                        (to sample from) with the following columns:
                                        'epitope'
                                        'binding'
        num_peptides                :   Number of peptides.
        num_peptides_immunogenic    :   Number of immunogenic peptides.
        peptide_length              :   Peptide length. Please note that
                                        all peptides sampled will be of this length.
        sampling_method             :   Sampling method. Allowed options:
                                        'random',
                                        'levenshtein',
                                        'alanine_scanning'

        Returns
        ------
        df_peptides                 :   Pandas DataFrame of sampled peptides
                                        with the following columns:
                                        'peptide_id'
                                        'epitope'
                                        'binding'
        """
        # Step to remove the HLA column from the reference and drop duplicates
        df_ref_peptides = df_ref_peptides.loc[:,['epitope', 'binding']].drop_duplicates()
        num_positives = num_peptides_immunogenic
        num_negatives = num_peptides - num_positives
        random_seed = Experiment.generate_random_seed()
        if sampling_method == 'random':
            # Sample the peptides from the reference dataframe without replacement
            df_positives = df_ref_peptides[df_ref_peptides['binding'] == 1].sample(n=num_positives, random_state=random_seed)
            df_negatives = df_ref_peptides[df_ref_peptides['binding'] == 0].sample(n=num_negatives, random_state=random_seed)
        elif sampling_method == 'levenshtein':
            # This is done similar to a sliding window based search for different neo-epitopes (LENS)
            # Pop a random row from the reference dataframe
            sampled_peptide = df_ref_peptides.sample(n=1, random_state=random_seed)
            df_ref_peptides.drop(sampled_peptide.index, inplace=True)
            # Sort peptides on distance to the sampled peptide from least to greatest Levenshtein distance
            df_ref_peptides['distance'] = df_ref_peptides['epitope'].apply(lambda x: Experiment.levenshtein_distance(x, sampled_peptide['epitope'].values[0]))
            df_ref_peptides.sort_values(by='distance', inplace=True)
            # Select the top n peptides from the sorted dataframe stratified by immunogenicity
            df_positives = df_ref_peptides[df_ref_peptides['binding'] == 1].head(num_positives)
            df_negatives = df_ref_peptides[df_ref_peptides['binding'] == 0].head(num_negatives)
        elif sampling_method == 'alanine_scanning':
            length_filter_mask = (df_ref_peptides['epitope'].str.len() == peptide_length)
            df_ref_peptide_length_filtered = df_ref_peptides.loc[length_filter_mask]

            # Sample positive peptides
            num_positives_ = num_positives
            sampled_peptides = set()
            positives_data = {
                'epitope': [],
                'binding': [],
                'reference_sequence': []
            }
            terminate = False
            while not terminate:
                df_positives = df_ref_peptide_length_filtered.loc[
                    (df_ref_peptide_length_filtered['binding'] == 1) &
                    ~(df_ref_peptide_length_filtered['epitope'].isin(sampled_peptides)),:
                ].sample(n=1, random_state=random_seed)
                if len(df_positives) == 0:
                    raise Exception('No more positive reference peptides to sample from')
                peptide_sequence = df_positives['epitope'].values.tolist()[0]
                sampled_peptides.add(peptide_sequence)
                peptide_sequences = Experiment.perform_alanine_scanning(peptide_sequence=peptide_sequence)
                if num_positives_ >= len(peptide_sequences):
                    for peptide_sequence_ in peptide_sequences:
                        positives_data['epitope'].append(peptide_sequence_)
                        positives_data['reference_sequence'].append(peptide_sequence)
                        # Logic for whether to flip the labels for the peptide based on the anchor residue
                        if peptide_sequence_[1] == 'A' and peptide_sequence[1] != 'A':
                            positives_data['binding'].append(0)
                            num_negatives -= 1
                        elif peptide_sequence_[8] == 'A' and peptide_sequence[8] != 'A':
                            positives_data['binding'].append(0)
                            num_negatives -= 1
                        else:
                            positives_data['binding'].append(1)
                            num_positives_ -= 1
                        if num_positives_ == 0:
                            terminate = True
                            break
                else:
                    # Randomize the remaining peptides to take such that not just the first k positions are alanine
                    positive_indices = random.sample(list(range(0, len(peptide_sequences))), num_positives_)
                    for idx in range(0, len(peptide_sequences)):
                        peptide_sequence_ = peptide_sequences[idx]
                        positives_data['epitope'].append(peptide_sequence_)
                        positives_data['reference_sequence'].append(peptide_sequence)
                        if idx in positive_indices:
                            if peptide_sequence_[1] == 'A' and peptide_sequence[1] != 'A':
                                positives_data['binding'].append(0)
                                num_negatives -= 1
                            elif peptide_sequence_[8] == 'A' and peptide_sequence[8] != 'A':
                                positives_data['binding'].append(0)
                                num_negatives -= 1
                            else:
                                positives_data['binding'].append(1)
                                num_positives_ -= 1
                        else:
                            positives_data['binding'].append(0)
                            num_negatives -= 1
                        if num_positives_ == 0:
                            terminate = True
                            break
            df_positives = pd.DataFrame(positives_data)

            # Sample negative peptides
            num_negatives_ = num_negatives
            sampled_peptides = set()
            negatives_data = {
                'epitope': [],
                'binding': [],
                'reference_sequence': []
            }
            terminate = False
            while not terminate:
                df_negatives = df_ref_peptide_length_filtered.loc[
                    (df_ref_peptide_length_filtered['binding'] == 0) &
                    ~(df_ref_peptide_length_filtered['epitope'].isin(sampled_peptides)),:
                ].sample(n=1, random_state=random_seed)
                if len(df_negatives) == 0:
                    raise Exception('No more negative reference peptides to sample from')
                peptide_sequence = df_negatives['epitope'].values.tolist()[0]
                sampled_peptides.add(peptide_sequence)
                peptide_sequences = Experiment.perform_alanine_scanning(peptide_sequence=peptide_sequence)
                for peptide_sequence_ in peptide_sequences:
                    negatives_data['epitope'].append(peptide_sequence_)
                    negatives_data['binding'].append(0)
                    negatives_data['reference_sequence'].append(peptide_sequence)
                    num_negatives_ -= 1
                    if num_negatives_ == 0:
                        terminate = True
                        break
            df_negatives = pd.DataFrame(negatives_data)
        else:
            raise Exception('Unsupported sampling_method: %s' % sampling_method)

        # Concatenate the positive and negative peptides and shuffle the dataframe
        df_peptides = pd.concat([df_positives, df_negatives], axis=0).sample(frac=1.0, random_state=random_seed)
        # Reset the index of the dataframe
        df_peptides = df_peptides.reset_index(drop=True)
        # Add a peptide_id column to the dataframe
        df_peptides['peptide_id'] = ['peptide_%i' % i for i in range(1, len(df_peptides) + 1)]
        if sampling_method == 'alanine_scanning':
            df_peptides = df_peptides.loc[:, ['peptide_id', 'epitope', 'binding', 'reference_sequence']]
        else:
            df_peptides = df_peptides.loc[:, ['peptide_id', 'epitope', 'binding']]
        return df_peptides

    @staticmethod
    def simulate_peptide_spot_counts(
            df_peptides: pd.DataFrame,
            num_coverage: int,
            random_effects: bool,
            mu_immunogenic: float,
            mu_nonimmunogenic: float,
            dispersion_factor: float
    ) -> pd.DataFrame:
        """
        Simulate peptide level spot counts.

        Parameters
        ----------
        df_peptides         :   Pandas DataFrame with the following columns:
                                'peptide_id'
                                'epitope'
                                'binding'
        num_coverage        :   Coverage number.
        random_effects      :   If True, enables random effects; uses
                                mu_immunogenic, mu_nonimmunogenic, and
                                dispersion_factor.
                                If False, mu_immunogenic = 1, mu_nonimmunogenic = 0,
                                dispersion_factor = 0.
        mu_immunogenic      :   Mean spot count for an immunogenic peptide.
        mu_nonimmunogenic   :   Mean spot count for a non-immunogenic peptide.
        dispersion_factor   :   Disperson factor.

        Returns
        -------
        df_peptides         :   Pandas DataFrame with the following columns:
                                'peptide_id'
                                'epitope'
                                'binding'
                                'coverage'
                                'spot_count'
        """
        # Step 1. Determine immunogenic and non-immunogenic mean and dispersion factor values
        if random_effects:
            # Simulate using random_effects (negative binomial distribution)
            immunogenic_mean = mu_immunogenic
            non_immunogenic_mean = mu_nonimmunogenic
        else:
            # Simulate using deterministic (binary) parameters
            immunogenic_mean = 1
            non_immunogenic_mean = 0
            dispersion_factor = 0

        # Step 2. Simulate peptide spot counts
        data = {
            'peptide_id': [],
            'epitope': [],
            'binding': [],
            'coverage': [],
            'spot_count': []
        }
        for index, row in df_peptides.iterrows():
            curr_id = row['peptide_id']
            curr_epitope = row['epitope']
            curr_binding = row['binding']
            if curr_binding == 1:
                mean_ = immunogenic_mean
            else:
                mean_ = non_immunogenic_mean
            peptide_spot_counts = Experiment.sample_spot_counts(
                mean=mean_,
                dispersion_factor=dispersion_factor,
                num_coverage=num_coverage
            )
            curr_coverage = 1
            for peptide_spot_count in peptide_spot_counts:
                data['peptide_id'].append(curr_id)
                data['epitope'].append(curr_epitope)
                data['binding'].append(curr_binding)
                data['coverage'].append(curr_coverage)
                data['spot_count'].append(peptide_spot_count)
                curr_coverage += 1
        return pd.DataFrame(data)

    @staticmethod
    def aggregate_pool_spot_counts(
            df_assignments: pd.DataFrame,
            df_peptides_spot_counts: pd.DataFrame,
            mu_nonimmunogenic: float,
            dispersion_factor: float,
            false_negative_rate: float,
            num_peptides_per_pool: int
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
        df_assignments          :  Pandas DataFrame with the following columns:
                                    'coverage_id'
                                    'pool_id'
                                    'peptide_id'
        df_peptides_spot_counts :   Pandas DataFrame with the following columns:
                                    'peptide_id'
                                    'epitope'
                                    'binding'
                                    'coverage'
                                    'spot_count'
        mu_nonimmunogenic       :   Non-immunogenic peptide spot count mean.
        dispersion_factor       :   Dispersion factor.
        false_negative_rate     :   Probability at which an immunogenic peptide
                                    spot count (at single coverage) becomes
                                    a non-immunogenic peptide spot count.
        num_peptides_per_pool   :   Number of peptides per pool.

        Returns
        -------
        df_readout              :   Pandas DataFrame with the following columns:
                                    'pool_id'
                                    'spot_count'
        """
        # Step 1. Aggregate peptide spot count to pool spot counts
        # key   = pool ID
        # value = pool spot count
        pool_spot_counts = {}
        for peptide_id in df_assignments['peptide_id'].unique():
            binding = df_peptides_spot_counts.loc[df_peptides_spot_counts['peptide_id'] == peptide_id, 'binding'].values.tolist()[0]
            pool_ids = df_assignments.loc[df_assignments['peptide_id'] == peptide_id, 'pool_id'].values.tolist()
            peptide_spot_counts = df_peptides_spot_counts.loc[df_peptides_spot_counts['peptide_id'] == peptide_id, 'spot_count'].values.tolist()
            random.shuffle(pool_ids)
            for i in range(0, len(pool_ids)):
                curr_pool_id = pool_ids[i]
                curr_peptide_spot_count = peptide_spot_counts[i]
                if curr_pool_id not in pool_spot_counts.keys():
                    pool_spot_counts[curr_pool_id] = 0
                if binding == 1:
                    # Randomize to include a false negative
                    if random.uniform(0.0, 1.0) <= false_negative_rate:
                        curr_peptide_spot_count = Experiment.sample_spot_counts(
                            mean=mu_nonimmunogenic,
                            dispersion_factor=dispersion_factor,
                            num_coverage=1
                        )[0]
                        curr_peptide_spot_count = float(curr_peptide_spot_count) / float(num_peptides_per_pool)
                else:
                    curr_peptide_spot_count = float(curr_peptide_spot_count) / float(num_peptides_per_pool)
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

    @staticmethod
    def identify_hit_pools(
            df_readouts: pd.DataFrame,
            num_peptides_per_pool: int,
            method: Literal['threshold', 'adaptive', 'combined'],
            threshold: float,
            mu_nonimmunogenic: float,
            dispersion_factor: float,
            alpha: float = 0.05,
            negative_control_sampling_iters: int = 10000,
    ) -> List[str]:
        """
        Identify the hit pools given the simulated spot counts for each pool.
        Pass the hit pools to the deconvolution method implemented by the
        different solvers. Based on the Empirical Rule (ER) response definition
        criteria from hte Moodie et al. (2010) paper: Response definition
        criteria for ELISPOT assays revisited:
        https://pubmed.ncbi.nlm.nih.gov/20549207/

        Parameters
        ----------
        df_readouts                         :   Pandas DataFrame with the following columns:
                                                'pool_id'
                                                'spot_count'
        num_peptides_per_pool               :   Number of peptides per pool in the ELISpot configuration.
        method                              :   Hit pool identification method.

                                                If 'threshold', hit pools are those with
                                                spot counts greater than or equal to the threshold
                                                (Dubey et al.).

                                                If 'adaptive', hit pools are identified based
                                                on the control pools (e.g. DMSO). Please note that
                                                'adaptive' only works whether random effects were
                                                applied to the peptide level spot counts.

                                                If 'combined', hit pools are identified based
                                                on the control pools (e.g. DMSO) and the threshold.
        threshold                           :   Threshold. If unspecified, the threshold is set
                                                by default based on whether random effects were
                                                applied to the peptide level spot counts.
        mu_nonimmunogenic                   :   Mean spot count for a non-immunogenic peptide.
        dispersion_factor                   :   Dispersion factor.
        alpha                               :   Significance level (default: 0.05).
        negative_control_sampling_iters     :   Number of times to sample negative control samples.

        Returns
        -------
        hit_pool_ids                        :   Hit pool IDs.
        """
        significance = (1 - alpha) * 100
        if method == 'threshold':
            hit_pool_ids = []
            for index, row in df_readouts.iterrows():
                pool_id = row['pool_id']
                spot_count = row['spot_count']
                if spot_count >= threshold:
                    hit_pool_ids.append(pool_id)
        elif method == 'adaptive':
            assert 0 < alpha < 1, "alpha must be between 0 and 1."
            # Simulate modified DMSO controls
            # Normal DMSO controls are simulated with the non-immunogenic mean across the whole pool
            # Here we simulate the DMSO controls with the non-immunogenic mean per peptide and add across the whole pool
            negative_controls = []
            for _ in range(0, negative_control_sampling_iters):
                negative_controls.append(sum(
                    Experiment.sample_spot_counts(
                        mean=mu_nonimmunogenic,
                        dispersion_factor=dispersion_factor,
                        num_coverage=num_peptides_per_pool
                )))
            hit_pool_ids = []
            for index, row in df_readouts.iterrows():
                pool_id = row['pool_id']
                spot_count = row['spot_count']
                if spot_count >= np.percentile(negative_controls, significance):
                    hit_pool_ids.append(pool_id)
        elif method == 'combined':
            thresholded_pool_ids = []
            for index, row in df_readouts.iterrows():
                pool_id = row['pool_id']
                spot_count = row['spot_count']
                if spot_count >= threshold:
                    thresholded_pool_ids.append(pool_id)
            negative_controls = []
            for _ in range(negative_control_sampling_iters):
                negative_controls.append(sum(
                    Experiment.sample_spot_counts(
                        mean=mu_nonimmunogenic,
                        dispersion_factor=dispersion_factor,
                        num_coverage=num_peptides_per_pool
                )))
            adaptive_pool_ids = []
            for index, row in df_readouts.iterrows():
                pool_id = row['pool_id']
                spot_count = row['spot_count']
                if spot_count >= np.percentile(negative_controls, significance):
                    adaptive_pool_ids.append(pool_id)
            hit_pool_ids = list(set(thresholded_pool_ids + adaptive_pool_ids))
        else:
            raise ValueError("Invalid method. Please choose from ['threshold', 'adaptive', or 'combined'].")
        return hit_pool_ids

    @staticmethod
    def evaluate_deconvolution_results(
            experiment_id: str,
            df_peptides: pd.DataFrame,
            block_assignment: BlockAssignment,
            deconvolution_result: DeconvolutionResult
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Evaluates deconvolution results.

        Parameters
        ----------
        experiment_id               :   Experiment ID.
        df_peptides                 :   Pandas DataFrame with the following columns (ground truth):
                                        'peptide_id'
                                        'epitope'
                                        'binding'
        block_assignment            :   BlockAssignment object.
        deconvolution_result        :   DeconvolutionResult object.

        Returns
        -------
        df_evaluation_metrics       :   Pandas DataFrame with the following columns:
                                        'experiment_id'
                                        'predicted_total_pools'
                                        'precision'
                                        'sensitivity'
                                        'specificity'
                                        'aucroc'
        df_evaluation_results       :   Pandas DataFrame with the following columns:
                                        'peptide_id'
                                        'binding'
                                        'peptide_spot_count'
        """
        df_assignments = block_assignment.to_dataframe()
        df_deconvolution_result = deconvolution_result.to_dataframe()
        confident_hit_peptide_ids = []
        candidate_hit_peptide_ids = []
        for _, row in df_deconvolution_result.iterrows():
            if row['deconvolution_result'] == DeconvolutionLabels.CANDIDATE_HIT:
                candidate_hit_peptide_ids.append(row['peptide_id'])
            if row['deconvolution_result'] == DeconvolutionLabels.CONFIDENT_HIT:
                confident_hit_peptide_ids.append(row['peptide_id'])

        # Step 1. Compute evaluation metrics
        predicted_total_pools = len(df_assignments['pool_id'].unique()) + len(candidate_hit_peptide_ids)
        hit_peptide_ids = confident_hit_peptide_ids + candidate_hit_peptide_ids
        no_hit_peptide_ids = list(set(df_assignments['peptide_id'].unique()) - set(hit_peptide_ids))
        positive_peptide_ids = df_peptides.loc[df_peptides['binding'] == 1, 'peptide_id'].unique()
        negative_peptide_ids = df_peptides.loc[df_peptides['binding'] == 0, 'peptide_id'].unique()
        if len(hit_peptide_ids) == 0:
            precision = 0
        else:
            precision = len(set(positive_peptide_ids).intersection(set(hit_peptide_ids))) / len(hit_peptide_ids)
        sensitivity = len(set(positive_peptide_ids).intersection(set(hit_peptide_ids))) / len(positive_peptide_ids)
        specificity = len(set(negative_peptide_ids).intersection(set(no_hit_peptide_ids))) / len(negative_peptide_ids)
        y_true = []
        y_pred = []
        results_data = {
            'peptide_id': [],
            'binding': [],
            'peptide_spot_count': []
        }
        for curr_peptide_id in df_peptides['peptide_id'].unique():
            curr_y_true = df_peptides.loc[df_peptides['peptide_id'] == curr_peptide_id, 'binding'].values.tolist()[0]
            df_matched = df_deconvolution_result.loc[df_deconvolution_result['peptide_id'] == curr_peptide_id,:]
            if len(df_matched) == 0:
                curr_y_pred = 0.0
            else:
                curr_y_pred = df_matched['estimated_peptide_spot_count'].values.tolist()[0]
            y_true.append(float(curr_y_true))
            y_pred.append(float(curr_y_pred))
            results_data['peptide_id'].append(curr_peptide_id)
            results_data['binding'].append(float(curr_y_true))
            results_data['peptide_spot_count'].append(float(curr_y_pred))
        df_evaluation_results = pd.DataFrame(results_data)
        aucroc = roc_auc_score(y_true, y_pred)
        metrics_data = {
            'experiment_id': [experiment_id],
            'predicted_total_pools': [predicted_total_pools],
            'precision': [precision],
            'sensitivity': [sensitivity],
            'specificity': [specificity],
            'aucroc': [aucroc]
        }
        df_evaluation_metrics = pd.DataFrame(metrics_data)
        return df_evaluation_metrics, df_evaluation_results

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
        # Iterate over the actual peptide sequence
        for i in range(0, len(peptide_sequence)):
            peptide_sequence_ = list(peptide_sequence)
            peptide_sequence_[i] = 'A'
            peptide_sequences.append(''.join(peptide_sequence_))
        return peptide_sequences

    @staticmethod
    def sample_spot_counts(
            mean: float,
            dispersion_factor: float,
            num_coverage: int
    ) -> List[float]:
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
        num_coverage        :   Coverage number.

        Returns
        -------
        spot_counts         :   List of integers. Each integer is a peptide level
                                spot count for a coverage.
        """
        if dispersion_factor == 0:
            return [mean] * num_coverage
        elif dispersion_factor < 1:
            raise ValueError("Dispersion_factor must be an integer greater than or equal to 1, or 0.")
        elif dispersion_factor == 1:
            return list(np.random.poisson(mean, num_coverage))
        else:
            variance = mean * dispersion_factor
            p = mean / variance
            r = mean ** 2 / (variance - mean)
            return list(np.random.negative_binomial(r, p, num_coverage))
