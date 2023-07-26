import pandas as pd
import os
import random
from acesim import Experiment, AceSolver, RandomizedAssignmentSolver


NUM_PROCESSES = 24
NUM_ITERATIONS = 100
REFERENCE_PEPTIDES_CSV_FILE = '/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/raw/iedb_mmer_all.csv'
DESIGNS_CSV_FILE = '/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/raw/ace_vs_naive_approach_alanine_scanning_designs.csv'
TRAINED_MODEL_FILE = '/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/raw/seq_sim_trained_model.pt'
OUTPUT_DIR = '/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/01_benchmark_ace_vs_naive_approach/alanine_scanning'
GOLFY_INIT_MODE = 'greedy'
GOLFY_MAX_ITERS = 2000


def get_solvers():
    random_seed = Experiment.generate_random_seed()
    ace_solver_1 = AceSolver(
        name='ace_golfy_clusteroff_noextrapools',
        cluster_peptides=False,
        random_seed=random_seed,
        mode='golfy',
        trained_model_file='',
        golfy_max_iters=GOLFY_MAX_ITERS,
        golfy_init_mode=GOLFY_INIT_MODE,
        golfy_allow_extra_pools=False
    )
    ace_solver_2 = AceSolver(
        name='ace_golfy_clusteron_noextrapools',
        cluster_peptides=True,
        random_seed=random_seed,
        mode='golfy',
        trained_model_file=TRAINED_MODEL_FILE,
        sim_threshold=0.7,
        sim_fxn='euclidean',
        golfy_max_iters=GOLFY_MAX_ITERS,
        golfy_init_mode=GOLFY_INIT_MODE,
        golfy_allow_extra_pools=False
    )
    random_solver = RandomizedAssignmentSolver(
        name='randomized_block_assignment',
        random_seed=random_seed
    )
    return [
        ace_solver_1,
        ace_solver_2,
        random_solver
    ]


if __name__ == "__main__":
    df_experiments = pd.read_csv(DESIGNS_CSV_FILE)
    df_results_all = pd.DataFrame()
    print('Started benchmarking ACE (golfy) vs naive approach')
    for index, row in df_experiments.iterrows():
        num_peptides = row['num_peptides']
        num_peptides_per_pool = row['num_peptides_per_pool']
        num_coverage = row['num_coverage']
        num_true_positive_peptides = row['num_true_positive_peptides']
        print('Started running experiment for the following configuration:')
        print('\tNumber of peptides: %i' % num_peptides)
        print("\tNumber of peptides per pool: %i" % num_peptides_per_pool)
        print("\tNumber of coverage: %i" % num_coverage)
        print("\tNumber of true positive peptides: %i" % num_true_positive_peptides)
        experiment_id = random.randint(1, 1000000000)
        experiment = Experiment(
            experiment_id=experiment_id,
            num_peptides=num_peptides,
            num_peptides_per_pool=num_peptides_per_pool,
            coverage=num_coverage,
            num_positives=num_true_positive_peptides,
            peptide_sampling_method='alanine_scanning',
            solvers=get_solvers(),
            df_ref_peptides=pd.read_csv(REFERENCE_PEPTIDES_CSV_FILE),
            num_processes=NUM_PROCESSES
        )
        df_results, df_assignments = experiment.run(num_iterations=NUM_ITERATIONS)
        output_dir = OUTPUT_DIR + '/%ipeptides_%iperpool_%ix' % (num_peptides, num_peptides_per_pool, num_coverage)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        df_results.drop(columns=['peptides'], inplace=True)
        df_results.to_csv(
            output_dir + '/ace_vs_naive_approach_benchmark_experiment_results_%i.tsv' % experiment_id,
            sep='\t',
            index=False
        )
        df_assignments.to_csv(
            output_dir + '/ace_vs_naive_approach_benchmark_experiment_assignments_%i.tsv' % experiment_id,
            sep='\t',
            index=False
        )
        print('Finished running the above experiment.')
    print('Finished benchmarking ACE (golfy) vs naive approach')
