import numpy as np
import pandas as pd
from golfy.evaluation import EvaluationResult
from golfy.optimization import improve_solution
from acelib.main import run_ace_deconvolve
from acelib.deconvolution import DeconvolutionLabels
from acesim.experiment import Experiment
from acesim.solver_repeated_design import RepeatedDesignSolver
from acesim.utilities import convert_block_assignment_to_golfy_design, convert_golfy_design_to_block_assignment


NUM_PEPTIDES = 400
MAX_PEPTIDES_PER_POOL = 20
COVERAGE = 3
NUM_POSITIVES = 4
SIM_ITERS = 100
OPTIM_ITERS = 100
REFERENCE_PEPTIDES_CSV_FILE = '/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/raw/held_out_data_w_negatives.csv'
OUTPUT_DIR = '/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/01_benchmark_ace/01_precision_recall_violations_iterations'


def optimize_block_assignments(
        n_peptides, 
        n_positives, 
        n_peptides_per_pool, 
        n_replicates, 
        n_optimization_iters
):
    # Step 1. Generate dummy peptides
    peptides = []
    for i in range(1, n_peptides + 1):
        peptides.append(('peptide_%i' % i, ''))

    # Step 2. Sample positive peptide IDs
    positive_peptide_ids = np.random.choice([peptide_id for peptide_id, _ in peptides], n_positives, replace=False)

    # Step 3. Sample the spot counts
    peptide_spot_counts = {}
    for peptide_id, _ in peptides:
        if peptide_id in positive_peptide_ids:
            peptide_spot_counts[peptide_id] = Experiment.sample_spot_counts(100, 1, n_replicates)
        else:
            peptide_spot_counts[peptide_id] = Experiment.sample_spot_counts(10, 1, n_replicates)

    # Step 4. Generate a repeated block assignment
    repeated_design_solver = RepeatedDesignSolver(name='')
    repeated_block_assignment, _ = repeated_design_solver.generate_assignment(
        peptides=peptides,
        num_peptides_per_pool=n_peptides_per_pool,
        num_coverage=n_replicates,
        random_seed=1
    )
    s, _ = convert_block_assignment_to_golfy_design(
        block_assignment=repeated_block_assignment
    )

    iter_block_designs = [repeated_block_assignment]
    for i in range(n_optimization_iters):
        improve_solution(s)
        iter_block_designs.append(convert_golfy_design_to_block_assignment(s))

    # Step 5. Evaluate the block designs per iteration
    history = []
    for i in range(0, len(iter_block_designs)):
        # Index into the current block assignment and solver
        curr_block_assignment = iter_block_designs[i]
        df_assignment = iter_block_designs[i].to_dataframe()
            
        # Compute the pool spot counts per pool assigned by the solver
        df_readout = Experiment.aggregate_pool_spot_counts(
            df_assignment=df_assignment,
            peptide_spot_counts=peptide_spot_counts
        )

        # Deconvolve the hit peptides
        deconvolution_result = run_ace_deconvolve(
            df_readout=df_readout,
            block_assignment=curr_block_assignment,
            mode='empirical',
            statistical_min_peptide_activity=10.0,
            empirical_min_coverage=n_replicates,
            empirical_min_spot_count=100 + 10*max(0, n_peptides_per_pool - 2),
            verbose=False
        ) 
        df_deconvolution_result = deconvolution_result.to_dataframe()           

        # Filter empirical deconvolution for self.coverage peptide activity level
        df_deconvolution_result = df_deconvolution_result.loc[
            df_deconvolution_result['peptide_activity_level'] == n_replicates,:
        ]
        candidate_peptide_ids = df_deconvolution_result.loc[df_deconvolution_result['deconvolution_result'] == DeconvolutionLabels.CANDIDATE_HIT, 'peptide_id'].values.tolist()
        num_total_pools_empirical = len(df_assignment['pool_id'].unique()) + len(candidate_peptide_ids)

        hit_peptide_ids = df_deconvolution_result.loc[
            (df_deconvolution_result['deconvolution_result'].isin([DeconvolutionLabels.CONFIDENT_HIT, DeconvolutionLabels.CANDIDATE_HIT])),
            'peptide_id'
        ].values.tolist()

        # Compute the precision, recall, and f1
        tp = len(set(hit_peptide_ids).intersection(set(positive_peptide_ids)))
        fp = len({p for p in hit_peptide_ids if p not in positive_peptide_ids})
        fn = len({p for p in positive_peptide_ids if p not in hit_peptide_ids})

        precision = tp / (tp + fp) if tp + fp > 0 else 0
        recall = tp / (tp + fn) if tp + fn > 0 else 0
        f1 = (
            2 * precision * recall / (precision + recall)
            if precision + recall > 0
            else 0
        )

        evr = EvaluationResult(
                precision=precision,
                recall=recall,
                f1=f1,
                num_pools=num_total_pools_empirical,
                num_violations=curr_block_assignment.num_violations()
            )
        history.append(evr)
    return history


if __name__ == "__main__":
    data = {
        'simulation_number':[],
        'iteration': [],
        'precision': [],
        'recall': [],
        'f1': [],
        'num_violations': []
    }
    for i in range(SIM_ITERS):
        history = optimize_block_assignments(
            n_peptides=NUM_PEPTIDES, 
            n_positives=NUM_POSITIVES, 
            n_peptides_per_pool=MAX_PEPTIDES_PER_POOL, 
            n_replicates=COVERAGE, 
            n_optimization_iters=OPTIM_ITERS
        )
        for iteration, evaluation in enumerate(history):
            data['simulation_number'].append(i)
            data['iteration'].append(iteration)
            data['precision'].append(evaluation.precision)
            data['recall'].append(evaluation.recall)
            data['f1'].append(evaluation.f1)
            data['num_violations'].append(evaluation.num_violations)
    df = pd.DataFrame(data)
    df.to_csv(OUTPUT_DIR + '/precision_recall_violations.csv', index=False)
