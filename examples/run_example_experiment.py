import pandas as pd
from acesim import Experiment, AceSolver, RandomizedDesignSolver, RepeatedDesignSolver


NUM_PROCESSES = 6
NUM_ITERATIONS = 100
REFERENCE_PEPTIDES_CSV_FILE = '/Users/leework/Documents/Research/projects/project_ace/data/raw/held_out_data_w_negatives.csv'
TRAINED_MODEL_FILE = '/Users/leework/Documents/Research/projects/project_ace/data/raw/trained_model_w_negatives.pt'


if __name__ == "__main__":
    ace_solver_1 = AceSolver(
        name='ace_golfy_extrapools',
        cluster_peptides=False,
        mode='golfy',
        golfy_allow_extra_pools=True,
        trained_model_file=''
    )
    ace_solver_2 = AceSolver(
        name='ace_golfy_noextrapools',
        cluster_peptides=False,
        mode='golfy',
        golfy_allow_extra_pools=False,
        trained_model_file=''
    )
    ace_solver_3 = AceSolver(
        name='ace_golfy_extrapools',
        cluster_peptides=True,
        mode='golfy',
        golfy_allow_extra_pools=False,
        trained_model_file=TRAINED_MODEL_FILE
    )
    randomized_solver = RandomizedDesignSolver(name='random')
    repeated_solver = RepeatedDesignSolver(name='repeated')
    experiment = Experiment(
        experiment_id=1,
        num_peptides=25,
        num_positives=5,
        num_peptides_per_pool=5,
        coverage=3,
        solvers=[ace_solver_1, ace_solver_2, ace_solver_3, randomized_solver, repeated_solver],
        random_effects=True,
        df_ref_peptides=pd.read_csv(REFERENCE_PEPTIDES_CSV_FILE),
        num_processes=NUM_PROCESSES
    )
    df_results, df_assignments = experiment.run(num_iterations=NUM_ITERATIONS)
    df_results.to_csv(
        'outputs/example_experiment_results.tsv',
        sep='\t',
        index=False
    )
    df_assignments.to_csv(
        'outputs/example_experiment_assignments.tsv',
        sep='\t',
        index=False
    )

