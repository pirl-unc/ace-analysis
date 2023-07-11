import pandas as pd
import pkg_resources
from acesim import Experiment, AceSolver, BogeySolver


NUM_PROCESSES = 94
NUM_ITERATIONS = 1000
REFERENCE_PEPTIDES_CSV_FILE = '/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/raw/iedb_mmer_all.csv'
CONFIGURATIONS_CSV_FILE = '/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/raw/ace_vs_naive_approach_configurations.csv'
OUTPUT_DIR = '/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/01_benchmark_ace_vs_naive_approach'


def get_solvers():
    random_seed = Experiment.generate_random_seed()
    ace_solver_1 = AceSolver(
        name='ace_1',
        cluster_peptides=True,
        random_seed=random_seed,
        mode='golfy',
        trained_model_file=pkg_resources.resource_filename('acelib', 'resources/models/seq_sim_trained_model.pt'),
        sim_threshold=0.7,
        sim_fxn='euclidean',
        golfy_max_iters=2000,
        golfy_init_mode='greedy'
    )
    ace_solver_2 = AceSolver(
        name='ace_2',
        cluster_peptides=False,
        random_seed=random_seed,
        mode='golfy',
        trained_model_file='',
        golfy_max_iters=2000,
        golfy_init_mode='greedy'
    )
    bogey_solver = BogeySolver(
        name='bogey',
        random_seed=random_seed
    )
    return [
        ace_solver_1,
        ace_solver_2,
        bogey_solver
    ]

if __name__ == "__main__":
    df_experiments = pd.read_csv(CONFIGURATIONS_CSV_FILE)
    df_experiments = df_experiments.loc[df_experiments['num_peptides'] == 120,:]
    df_results_all = pd.DataFrame()
    print('Started benchmarking ACE vs naive approach')
    for index, row in df_experiments.iterrows():
        num_peptides = row['num_peptides']
        num_peptides_per_pool = row['num_peptides_per_pool']
        num_coverage = row['num_coverage']
        num_true_positive_peptides = row['num_true_positive_peptides']
        print('Started running experiment for the following configuration:')
        print('Number of peptides: %i' % num_peptides)
        print("Number of peptides per pool: %i" % num_peptides_per_pool)
        print("Number of coverage: %i" % num_coverage)
        print("Number of true positive peptides: %i" % num_true_positive_peptides)
        experiment = Experiment(
            num_peptides=num_peptides,
            num_peptides_per_pool=num_peptides_per_pool,
            coverage=num_coverage,
            peptide_scan=True,
            num_positives=num_true_positive_peptides,
            solvers=get_solvers(),
            df_ref_peptides=pd.read_csv(REFERENCE_PEPTIDES_CSV_FILE),
            num_processes=NUM_PROCESSES
        )
        df_results = experiment.run(num_iterations=NUM_ITERATIONS)
        df_results['num_peptides'] = num_peptides
        df_results['num_peptides_per_pool'] = num_peptides_per_pool
        df_results['num_coverage'] = num_coverage
        df_results['num_true_positive_peptides'] = num_true_positive_peptides
        df_results_all = pd.concat([df_results_all, df_results])
        print('Finished running the above experiment.')
    df_results_all.to_csv(OUTPUT_DIR + '/ace_vs_naive_approach_benchmark_experiment_results.csv', index=False)
    print('Finished benchmarking ACE vs naive approach')
