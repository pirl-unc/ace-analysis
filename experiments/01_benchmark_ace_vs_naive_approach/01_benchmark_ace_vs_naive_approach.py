import pandas as pd
import pkg_resources
from acesim import Experiment, AceSolver, BogeySolver


NUM_PROCESSES = 4
NUM_ITERATIONS = 100
REFERENCE_PEPTIDES_CSV_FILE = '/Users/leework/Documents/Research/projects/project_ace/scripts/ace-analysis/test/data/iedb_mmer_all.csv'
OUTPUT_DIR = '/Users/leework/Documents/Research/projects/project_ace/data/processed/03_benchmark_ace_vs_naive_approach'


if __name__ == "__main__":
    ace_solver = AceSolver(
        name='ace',
        cluster_peptides=True,
        random_seed=1,
        mode='golfy',
        trained_model_file=pkg_resources.resource_filename('acelib', 'resources/models/seq_sim_trained_model.pt')
    )
    bogey_solver = BogeySolver(
        name='bogey'
    )
    experiment = Experiment(
        num_peptides=120,
        num_peptides_per_pool=12,
        coverage=3,
        num_positives=10,
        solvers=[ace_solver, bogey_solver],
        df_ref_peptides=pd.read_csv(REFERENCE_PEPTIDES_CSV_FILE),
        num_processes=NUM_PROCESSES
    )
    df_results = experiment.run(num_iterations=NUM_ITERATIONS)
    df_results.to_csv(OUTPUT_DIR + '/ace_vs_naive_approach_benchmark_experiment_results.csv', index=False)


