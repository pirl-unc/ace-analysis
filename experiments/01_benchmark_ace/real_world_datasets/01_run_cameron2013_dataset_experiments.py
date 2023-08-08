import pandas as pd
import random
import os
from acesim import Experiment, AceSolver, RandomizedDesignSolver, RepeatedDesignSolver


NUM_PROCESSES = 64
NUM_ITERATIONS = 100
REFERENCE_PEPTIDES_CSV_FILE = '/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/raw/cameron_et_al_sci_trans_med_2013.csv'
TRAINED_MODEL_FILE = '/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/raw/seq_sim_trained_model.pt'
OUTPUT_DIR = '/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/01_benchmark_ace/real_world_datasets/cameron_et_al_sci_trans_med_2013'
GOLFY_INIT_MODE = 'greedy'
GOLFY_MAX_ITERS = 2000
MIN_PEPTIDE_ACTIVITY = 10.0
RANDOM_EFFECTS = True


def get_solvers():
    ace_solver_1 = AceSolver(
        name='ace_golfy_clusteroff_noextrapools',
        cluster_peptides=False,
        mode='golfy',
        trained_model_file='',
        golfy_max_iters=GOLFY_MAX_ITERS,
        golfy_init_mode=GOLFY_INIT_MODE,
        golfy_allow_extra_pools=False
    )
    ace_solver_2 = AceSolver(
        name='ace_golfy_clusteron_noextrapools',
        cluster_peptides=True,
        mode='golfy',
        trained_model_file=TRAINED_MODEL_FILE,
        sim_threshold=0.8,
        sim_fxn='euclidean',
        golfy_max_iters=GOLFY_MAX_ITERS,
        golfy_init_mode=GOLFY_INIT_MODE,
        golfy_allow_extra_pools=False
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
    print('Started running an experiment on Cameron et al., Sci Trans Med 2013, dataset')
    num_peptides = 36
    num_peptides_per_pool = 6
    num_coverage = 3
    num_true_positive_peptides = 3
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
        peptide_sampling_method='',
        solvers=get_solvers(),
        min_peptide_activity=MIN_PEPTIDE_ACTIVITY,
        dispersion_factor=4.0,
        random_effects=RANDOM_EFFECTS,
        df_ref_peptides=pd.read_csv(REFERENCE_PEPTIDES_CSV_FILE),
        num_processes=NUM_PROCESSES
    )
    df_results, df_assignments = experiment.run(num_iterations=NUM_ITERATIONS)
    output_dir = OUTPUT_DIR + '/%ipeptides_%iperpool_%ix' % (num_peptides, num_peptides_per_pool, num_coverage)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    df_results.drop(columns=['peptides'], inplace=True)
    df_results.to_csv(
        output_dir + '/cameron_et_al_sci_trans_med_2013_experiment_results_%i.tsv' % experiment_id,
        sep='\t',
        index=False
    )
    df_assignments.to_csv(
        output_dir + '/cameron_et_al_sci_trans_med_2013_experiment_assignments_%i.tsv' % experiment_id,
        sep='\t',
        index=False
    )
    print('Finished running the above experiment.')
    print('Finished running an experiment on Cameron et al., Sci Trans Med 2013, dataset')