import pandas as pd
import random
from acesim import Experiment, AceSolver, RandomizedDesignSolver, RepeatedDesignSolver


NUM_PROCESSES = 6
NUM_ITERATIONS = 100
REFERENCE_PEPTIDES_CSV_FILE = '/Users/leework/Documents/Research/projects/project_ace/data/raw/held_out_data_w_negatives.csv'
TRAINED_MODEL_FILE = '/Users/leework/Documents/Research/projects/project_ace/data/raw/trained_model_w_negatives.pt'


if __name__ == "__main__":
    ace_solver_1 = AceSolver(
        name='ace_golfy_clusteron',
        cluster_peptides=True,
        mode='golfy',
        golfy_allow_extra_pools=False,
        trained_model_file=TRAINED_MODEL_FILE
    )
    ace_solver_2 = AceSolver(
        name='ace_golfy_clusteroff',
        cluster_peptides=False,
        mode='golfy',
        golfy_allow_extra_pools=False,
        trained_model_file=''
    )
    randomized_solver = RandomizedDesignSolver(name='random')
    repeated_solver = RepeatedDesignSolver(name='repeated')
    experiment = Experiment(
        experiment_id=1,
        num_peptides=90,
        num_peptides_per_pool=9,
        coverage=3,
        num_positives=9,
        peptide_sampling_method='alanine_scanning',
        random_effects=True,
        solvers=[ace_solver_1, ace_solver_2, randomized_solver, repeated_solver],
        df_ref_peptides=pd.read_csv(REFERENCE_PEPTIDES_CSV_FILE),
        num_processes=NUM_PROCESSES
    )
    df_results, df_assignments = experiment.run(num_iterations=NUM_ITERATIONS)
    df_results.to_csv(
        'outputs/example_alanine_scanning_results.tsv',
        sep='\t',
        index=False
    )
    df_assignments.to_csv(
        'outputs/example_alanine_scanning_assignments.tsv',
        sep='\t',
        index=False
    )


