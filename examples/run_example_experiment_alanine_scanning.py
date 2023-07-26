import pandas as pd
import random
from acesim import Experiment, AceSolver, RandomizedDesignSolver, RepeatedDesignSolver


NUM_PROCESSES = 6
NUM_ITERATIONS = 10
REFERENCE_PEPTIDES_CSV_FILE = '/Users/leework/Documents/Research/projects/project_ace/scripts/ace-analysis/test/data/iedb_mmer_all.csv'


if __name__ == "__main__":
    ace_solver_1 = AceSolver(
        name='ace_golfy_clusteron',
        cluster_peptides=True,
        random_seed=Experiment.generate_random_seed(),
        mode='golfy',
        golfy_allow_extra_pools=False,
        trained_model_file='/Users/leework/Documents/Research/projects/project_ace/scripts/ace-analysis/notebooks/models/seq_sim_trained_model.pt'
    )
    ace_solver_2 = AceSolver(
        name='ace_golfy_clusteroff',
        cluster_peptides=False,
        random_seed=Experiment.generate_random_seed(),
        mode='golfy',
        golfy_allow_extra_pools=False,
        trained_model_file='/Users/leework/Documents/Research/projects/project_ace/scripts/ace-analysis/notebooks/models/seq_sim_trained_model.pt'
    )
    randomized_solver = RandomizedDesignSolver(name='random', random_seed=1)
    repeated_solver = RepeatedDesignSolver(name='repeated')
    experiment = Experiment(
        experiment_id=random.randint(1,1000000),
        num_peptides=90,
        num_peptides_per_pool=9,
        coverage=3,
        num_positives=9,
        peptide_sampling_method='alanine_scanning',
        solvers=[ace_solver_1, ace_solver_2, randomized_solver, repeated_solver],
        df_ref_peptides=pd.read_csv(REFERENCE_PEPTIDES_CSV_FILE),
        num_processes=NUM_PROCESSES
    )
    df_results, df_assignments = experiment.run(num_iterations=NUM_ITERATIONS)
    df_results.to_csv(
        'alanine_scanning_results.tsv',
        sep='\t',
        index=False
    )
    df_assignments.to_csv(
        'alanine_scanning_assignments.tsv',
        sep='\t',
        index=False
    )


