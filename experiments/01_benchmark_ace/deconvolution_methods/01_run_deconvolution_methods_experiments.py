import pandas as pd
import random
from acesim import Experiment, PrecomputedSolver
from acelib.block_assignment import BlockAssignment


NUM_PROCESSES = 64
NUM_ITERATIONS = 1000
REFERENCE_PEPTIDES_CSV_FILE = '/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/raw/iedb_mmer_all.csv'
RESOURCE_DIR = '/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/04_generate_resource_designs'
RESOURCE_CSV_FILE = '/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/raw/deconvolution_methods_resources.csv'
OUTPUT_DIR = '/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/01_benchmark_ace/deconvolution_methods'
MIN_PEPTIDE_ACTIVITY = 10.0
RANDOM_EFFECTS = True


if __name__ == "__main__":
    df_results_all = pd.DataFrame()
    df_resources = pd.read_csv(RESOURCE_CSV_FILE)
    experiment_id = 1
    for index, row in df_resources.iterrows():
        excel_file = RESOURCE_DIR + '/' + row['design_excel_file']
        num_peptides = row['num_peptides']
        num_peptides_per_pool = row['num_peptides_per_pool']
        num_coverage = row['num_coverage']
        num_positives = row['num_positives']
        block_assignment = BlockAssignment.read_excel_file(
            excel_file=excel_file,
            sheet_name='block_assignment'
        )
        precomputed_solver = PrecomputedSolver(
            name='precomputed',
            block_assignment=block_assignment
        )
        experiment = Experiment(
            experiment_id=experiment_id,
            num_peptides=num_peptides,
            num_peptides_per_pool=num_peptides_per_pool,
            coverage=num_coverage,
            num_positives=num_positives,
            solvers=[precomputed_solver],
            min_peptide_activity=MIN_PEPTIDE_ACTIVITY,
            random_effects=RANDOM_EFFECTS,
            df_ref_peptides=pd.read_csv(REFERENCE_PEPTIDES_CSV_FILE),
            num_processes=NUM_PROCESSES
        )
        df_results, df_assignments = experiment.run(num_iterations=NUM_ITERATIONS)
        df_results_all = pd.concat([df_results_all, df_results])
        experiment_id += 1

    df_results_all.to_csv(
        OUTPUT_DIR + '/deconvolution_methods_experiment_results.tsv',
        sep='\t',
        index=False
    )
