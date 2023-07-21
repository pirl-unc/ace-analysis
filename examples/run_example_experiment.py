import pandas as pd
from acesim import Experiment, AceSolver, BogeySolver


NUM_PROCESSES = 6
NUM_ITERATIONS = 100
REFERENCE_PEPTIDES_CSV_FILE = '/Users/leework/Documents/Research/projects/project_ace/scripts/ace-analysis/test/data/iedb_mmer_all.csv'
OUTPUT_DIR = '/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace_vs_naive_approach'


if __name__ == "__main__":
    ace_solver_1 = AceSolver(
        name='ace_golfy',
        cluster_peptides=True,
        random_seed=1,
        mode='golfy',
        trained_model_file='/Users/leework/Documents/Research/projects/project_ace/scripts/ace-analysis/notebooks/models/seq_sim_trained_model.pt'
    )
    ace_solver_2 = AceSolver(
        name='ace_sat',
        cluster_peptides=True,
        random_seed=1,
        mode='sat_solver',
        trained_model_file='/Users/leework/Documents/Research/projects/project_ace/scripts/ace-analysis/notebooks/models/seq_sim_trained_model.pt'
    )
    bogey_solver = BogeySolver(
        name='bogey',
        random_seed=1
    )
    experiment = Experiment(
        num_peptides=120,
        num_peptides_per_pool=12,
        coverage=3,
        num_positives=10,
        solvers=[ace_solver_1, ace_solver_2, bogey_solver],
        df_ref_peptides=pd.read_csv(REFERENCE_PEPTIDES_CSV_FILE),
        num_processes=NUM_PROCESSES
    )
    df_results = experiment.run(num_iterations=NUM_ITERATIONS)
    df_results.to_csv(OUTPUT_DIR + '/ace_vs_naive_approach_benchmark_experiment_results.csv', index=False)


