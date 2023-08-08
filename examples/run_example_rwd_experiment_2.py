import pandas as pd
from acesim import Experiment, AceSolver, RandomizedDesignSolver, RepeatedDesignSolver


NUM_PROCESSES = 6
NUM_ITERATIONS = 100
REFERENCE_PEPTIDES_CSV_FILE = '/Users/leework/Documents/Research/projects/project_ace/data/raw/tiling_window_covid_spike_9mers.csv'
TRAINED_MODEL_FILE = '/Users/leework/Documents/Research/projects/project_ace/data/raw/trained_model_w_negatives.pt'
GOLFY_INIT_MODE = 'greedy'
GOLFY_MAX_ITERS = 2000
MIN_PEPTIDE_ACTIVITY = 10.0
RANDOM_EFFECTS = True


def get_solvers():
    ace_solver_1 = AceSolver(
        name='ace_golfy_clusteroff_extrapools',
        cluster_peptides=False,
        mode='golfy',
        trained_model_file='',
        golfy_max_iters=GOLFY_MAX_ITERS,
        golfy_init_mode=GOLFY_INIT_MODE,
        golfy_allow_extra_pools=True
    )
    ace_solver_2 = AceSolver(
        name='ace_golfy_clusteron_extrapools',
        cluster_peptides=True,
        mode='golfy',
        trained_model_file=TRAINED_MODEL_FILE,
        sim_threshold=0.8,
        sim_fxn='euclidean',
        golfy_max_iters=GOLFY_MAX_ITERS,
        golfy_init_mode=GOLFY_INIT_MODE,
        golfy_allow_extra_pools=True
    )
    randomized_solver = RandomizedDesignSolver(name='randomized_block_assignment')
    repeated_solver = RepeatedDesignSolver(name='repeated_block_assignment')
    return [
        ace_solver_1,
        ace_solver_2,
        randomized_solver,
        repeated_solver
    ]


if __name__ == "__main__":
    df_results_all = pd.DataFrame()
    num_peptides = 1265
    num_peptides_per_pool = 35
    num_coverage = 3
    num_true_positive_peptides = 105
    experiment = Experiment(
        experiment_id=1,
        num_peptides=num_peptides,
        num_peptides_per_pool=num_peptides_per_pool,
        coverage=num_coverage,
        num_positives=num_true_positive_peptides,
        peptide_sampling_method='',
        solvers=get_solvers(),
        min_peptide_activity=MIN_PEPTIDE_ACTIVITY,
        random_effects=RANDOM_EFFECTS,
        df_ref_peptides=pd.read_csv(REFERENCE_PEPTIDES_CSV_FILE),
        num_processes=NUM_PROCESSES
    )
    df_results, df_assignments = experiment.run(num_iterations=NUM_ITERATIONS)
    df_results.drop(columns=['peptides'], inplace=True)
    df_results.to_csv(
        'outputs/covid_9mer_experiment_results.tsv',
        sep='\t',
        index=False
    )
    df_assignments.to_csv(
        'outputs/covid_9mer_experiment_assignments.tsv',
        sep='\t',
        index=False
    )
