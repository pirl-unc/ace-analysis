import pandas as pd
import random
import os
from acesim import Experiment, AceSolver, RandomizedDesignSolver, RepeatedDesignSolver


NUM_PROCESSES = 64
NUM_ITERATIONS = 1000
REFERENCE_PEPTIDES_CSV_FILE = '/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/raw/tan_et_al_infection_genetics_evolution_2021.csv'
TRAINED_MODEL_FILE = '/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/raw/trained_model_w_data_augmentation_b3000.pt'
OUTPUT_DIR = '/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/01_benchmark_ace/06_real_world_datasets/tan_et_al_inf_gen_evol_2021/noextrapools'
RANDOM_EFFECTS = True


def get_solvers():
    ace_solver_1 = AceSolver(
        name='ace_golfy_clusteroff_noextrapools',
        cluster_peptides=False,
        mode='golfy',
        trained_model_file='',
        golfy_allow_extra_pools=False
    )
    ace_solver_2 = AceSolver(
        name='ace_golfy_clusteron_noextrapools',
        cluster_peptides=True,
        mode='golfy',
        trained_model_file=TRAINED_MODEL_FILE,
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
    print('Started running an experiment on Tan et al. Infection Genetics Evolution 2021 dataset')
    num_peptides = 177
    num_peptides_per_pool = 10
    num_coverage = 3
    num_true_positive_peptides = 22
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
        '%s/tan_et_al_inf_gen_evol_2021_experiment_results_%i.tsv' % (output_dir, experiment_id),
        sep='\t',
        index=False
    )
    df_assignments.to_csv(
        '%s/tan_et_al_inf_gen_evol_2021_experiment_assignments_%i.tsv' % (output_dir, experiment_id),
        sep='\t',
        index=False
    )
    print('Finished running the above experiment.')
    print('Finished running an experiment on Tan et al. Infection Genetics Evolution 2021 dataset')
