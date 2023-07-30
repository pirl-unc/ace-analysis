import pandas as pd
import random
from acesim import Experiment, PrecomputedSolver, RandomizedDesignSolver, RepeatedDesignSolver
from acelib.block_assignment import BlockAssignment


NUM_PROCESSES = 6
NUM_ITERATIONS = 1000
REFERENCE_PEPTIDES_CSV_FILE = '/Users/leework/Documents/Research/projects/project_ace/scripts/ace-analysis/test/data/iedb_mmer_all.csv'
DESIGN_EXCEL_FILE = '/Users/leework/Documents/Research/projects/project_ace/data/processed/04_generate_resource_designs/100peptides_10perpool_3x.xlsx'


if __name__ == "__main__":
    block_assignment = BlockAssignment.read_excel_file(
        excel_file=DESIGN_EXCEL_FILE,
        sheet_name='block_assignment'
    )
    precomputed_solver = PrecomputedSolver(
        name='precomputed',
        block_assignment=block_assignment
    )
    experiment = Experiment(
        experiment_id=random.randint(1,1000000),
        num_peptides=100,
        num_peptides_per_pool=10,
        coverage=3,
        num_positives=10,
        solvers=[precomputed_solver],
        min_peptide_activity=10.0,
        random_effects=False,
        df_ref_peptides=pd.read_csv(REFERENCE_PEPTIDES_CSV_FILE),
        num_processes=NUM_PROCESSES
    )
    df_results, df_assignments = experiment.run(num_iterations=NUM_ITERATIONS)
    df_results.to_csv(
        'test_deconvolution_10percent_randomeffectsoff_results.tsv',
        sep='\t',
        index=False
    )
    df_assignments.to_csv(
        'test_deconvolution_10percent_randomeffectsoff_assignments.tsv',
        sep='\t',
        index=False
    )

