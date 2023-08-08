import pandas as pd
import random
from acesim import Experiment, PrecomputedSolver
from acelib.block_assignment import BlockAssignment


NUM_PROCESSES = 6
NUM_ITERATIONS = 100
REFERENCE_PEPTIDES_CSV_FILE = '/Users/leework/Documents/Research/projects/project_ace/data/raw/held_out_data_w_negatives.csv'
DESIGN_EXCEL_FILE = '/Users/leework/Documents/Research/projects/project_ace/scripts/ace-analysis/test/data/100peptides_10perpool_3x.xlsx'


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
        experiment_id=1,
        num_peptides=100,
        num_peptides_per_pool=10,
        coverage=3,
        num_positives=1,
        solvers=[precomputed_solver],
        random_effects=True,
        df_ref_peptides=pd.read_csv(REFERENCE_PEPTIDES_CSV_FILE),
        num_processes=NUM_PROCESSES
    )
    df_results, df_assignments = experiment.run(num_iterations=NUM_ITERATIONS)
    df_results.to_csv(
        'outputs/example_deconvolution_results.tsv',
        sep='\t',
        index=False
    )
    df_assignments.to_csv(
        'outputs/example_deconvolution_assignments.tsv',
        sep='\t',
        index=False
    )

