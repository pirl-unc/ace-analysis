import multiprocessing as mp
from golfy.optimization import improve_solution
from acesim import Experiment, RepeatedDesignConfigurator, AceDeconvolver
from acesim.utilities import *


NUM_PROCESSES = 64
NUM_PEPTIDES = 400
NUM_POSITIVES = 4
NUM_PEPTIDES_PER_POOL = 20
NUM_COVERAGE = 3
PEPTIDE_LENGTH = 9
MU_IMMUNOGENIC = 300.0
MU_NONIMMUNOGENIC = 5.0
DISPERSION_FACTOR = 0.0
FALSE_NEGATIVE_RATE = 0.0
PEPTIDE_SAMPLING_METHOD = 'random'
RANDOM_EFFECTS = True
DECONVOLUTION_METHOD = 'em'
SIM_ITERS = 100
OPTIM_ITERS = 100
REFERENCE_PEPTIDES_CSV_FILE = '/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/raw/references/held_out_data_w_negatives.csv'
OUTPUT_DIR = '/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/01_precision_recall_violations_iterations'


def optimize_block_assignments(
        experiment_id: str,
        block_assignment: BlockAssignment,
        df_peptides: pd.DataFrame,
        df_peptides_spot_counts: pd.DataFrame,
        mu_nonimmunogenic: float,
        dispersion_factor: float,
        false_negative_rate: float,
        deconvolution_method: str,
        min_pool_spot_count: float,
        num_coverage: int,
        num_optimization_iters: int
):
    """
    Optimizes a repeated block design iteratively.

    Parameters
    ----------
    df_peptides     :   Pandas DataFrame with the following columns:
                        'peptide_id',
                        'epitope',
                        'binding'

    """
    # Step 1. Convert the BlockAssignment to golfy design
    s, _ = convert_block_assignment_to_golfy_design(block_assignment=block_assignment)

    # Step 2. Improve the BlockAssignment
    iter_block_designs = [block_assignment]
    for _ in range(num_optimization_iters):
        improve_solution(s)
        iter_block_designs.append(convert_golfy_design_to_block_assignment(s))

    # Step 2. Evaluate the block designs per iteration
    df_evaluation_metrics_history = pd.DataFrame()
    for i in range(0, len(iter_block_designs)):
        # Index into the current block assignment and solver
        curr_block_assignment = iter_block_designs[i]
        df_assignments = iter_block_designs[i].to_dataframe()

        # Compute the pool spot counts per pool assigned by the solver
        df_readout = Experiment.aggregate_pool_spot_counts(
            df_assignments=df_assignments,
            df_peptides_spot_counts=df_peptides_spot_counts,
            mu_nonimmunogenic=mu_nonimmunogenic,
            dispersion_factor=dispersion_factor,
            false_negative_rate=false_negative_rate
        )

        # Deconvolve the hit peptides
        ace_deconvolver = AceDeconvolver(name='')
        deconvolution_result = ace_deconvolver.deconvolve(
            df_readout=df_readout,
            block_assignment=curr_block_assignment,
            method=deconvolution_method,
            min_coverage=num_coverage,
            min_pool_spot_count=min_pool_spot_count
        )

        # Evaluate
        df_evaluation_metrics, df_evaluation_results = Experiment.evaluate_deconvolution_results(
            experiment_id=experiment_id,
            df_peptides=df_peptides,
            block_assignment=curr_block_assignment,
            deconvolution_result=deconvolution_result
        )

        df_evaluation_metrics['iteration'] = [i]
        df_evaluation_metrics['num_violations'] = [curr_block_assignment.num_violations()]
        df_evaluation_metrics_history = pd.concat([df_evaluation_metrics_history, df_evaluation_metrics], axis=0, ignore_index=True)

    return df_evaluation_metrics_history


if __name__ == "__main__":
    # Step 1. Sample peptides
    df_peptides = Experiment.sample_peptides(
        df_ref_peptides=pd.read_csv(REFERENCE_PEPTIDES_CSV_FILE),
        num_peptides=NUM_PEPTIDES,
        num_peptides_immunogenic=NUM_POSITIVES,
        peptide_length=PEPTIDE_LENGTH,
        sampling_method=PEPTIDE_SAMPLING_METHOD
    )

    # Step 2. Simulate peptide spot counts
    df_peptides_spot_counts = Experiment.simulate_peptide_spot_counts(
        df_peptides=df_peptides,
        num_coverage=NUM_COVERAGE,
        random_effects=RANDOM_EFFECTS,
        mu_immunogenic=MU_IMMUNOGENIC,
        mu_nonimmunogenic=MU_NONIMMUNOGENIC,
        dispersion_factor=DISPERSION_FACTOR
    )

    df_peptides.to_csv(OUTPUT_DIR + "/sampled_peptides_ground_truth.tsv", sep='\t', index=False)
    df_peptides_spot_counts.to_csv(OUTPUT_DIR + "/sampled_peptides_spot_counts_ground_truth.tsv", sep='\t', index=False)

    # Step 3. Generate a repeated block design
    peptides = convert_peptides_dataframe_to_peptides(df_peptides=df_peptides)
    repeated_design_solver = RepeatedDesignConfigurator(name='repeated')
    repeated_block_assignment, _ = repeated_design_solver.generate_assignment(
        peptides=peptides,
        num_peptides_per_pool=NUM_PEPTIDES_PER_POOL,
        num_coverage=NUM_COVERAGE,
        plate_size=96
    )

    # Step 4. Improve a repeated block design
    pool = mp.Pool(processes=NUM_PROCESSES)
    results_async = []
    for i in range(0, SIM_ITERS):
        result_async = pool.apply_async(optimize_block_assignments, args=[
            'experiment_%i' % i,
            repeated_block_assignment,
            df_peptides,
            df_peptides_spot_counts,
            MU_NONIMMUNOGENIC,
            DISPERSION_FACTOR,
            FALSE_NEGATIVE_RATE,
            DECONVOLUTION_METHOD,
            MU_IMMUNOGENIC,
            NUM_COVERAGE,
            OPTIM_ITERS
        ])
        results_async.append(result_async)
    results = [result_async.get() for result_async in results_async]
    pool.close()
    pool.join()

    # Step 5. Write evaluation results
    df_evaluation_history = pd.DataFrame()
    for df_evaluation in results:
        df_evaluation_history = pd.concat([df_evaluation_history, df_evaluation], axis=0, ignore_index=True)
    df_evaluation_history.to_csv(OUTPUT_DIR + '/evaluations.tsv', sep='\t', index=False)
