#!/usr/bin/python3

"""
The purpose of this python3 script is to generate all ACE resource designs.

Last updated date: July 29, 2023
"""


import pandas as pd
from acelib.block_assignment import BlockAssignment
from acelib.block_design import BlockDesign
from acelib.constants import PlateTypes
from acelib.utilities import convert_dataframe_to_peptides
from acelib.main import  run_ace_sat_solver


DESIGN_CSV_FILE = '/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/raw/resource_designs.csv'
OUTPUT_DIR = '/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/04_generate_resource_designs'
NUM_PROCESSES = 32
SHUFFLE_ITERS = 1000
MAX_PEPTIDES_PER_BLOCK = 100
MAX_PEPTIDES_PER_POOL = 10


def generate_peptides_dataframe(num_peptides: int):
    data = {
        'peptide_id': ['peptide_%i' % i for i in range(1, num_peptides + 1)],
        'peptide_sequence': ['' for _ in range(0, num_peptides)]
    }
    return pd.DataFrame(data)


if __name__ == '__main__':
    df_designs = pd.read_csv(DESIGN_CSV_FILE)
    data = {
        'num_peptides': [],
        'num_peptides_per_pool': [],
        'num_coverage': [],
        'optimal': []
    }
    for index, value in df_designs.iterrows():
        num_peptides = value['num_peptides']
        num_peptides_per_pool = value['num_peptides_per_pool']
        num_coverage = value['num_coverage']
        print('Started generating: %i peptides / %i peptides per pool / %ix coverage' %
              (num_peptides, num_peptides_per_pool, num_coverage))
        df_peptides = generate_peptides_dataframe(num_peptides=num_peptides)
        block_design = BlockDesign(
            peptides=convert_dataframe_to_peptides(df_peptides=df_peptides),
            num_peptides_per_pool=num_peptides_per_pool,
            num_coverage=num_coverage,
            max_peptides_per_block=MAX_PEPTIDES_PER_BLOCK,
            disallowed_peptide_pairs=[],
            preferred_peptide_pairs=[]
        )
        block_assignment = run_ace_sat_solver(
            block_design=block_design,
            max_peptides_per_pool=MAX_PEPTIDES_PER_POOL,
            num_processes=NUM_PROCESSES,
            shuffle_iters=SHUFFLE_ITERS,
            verbose=True
        )
        df_assignment = block_assignment.to_dataframe()
        data['num_peptides'].append(num_peptides)
        data['num_peptides_per_pool'].append(num_peptides_per_pool)
        data['num_coverage'].append(num_coverage)
        data['optimal'].append(block_assignment.is_optimal(
            num_coverage=num_coverage,
            num_peptides_per_pool=num_peptides_per_pool,
            verbose=True
        ))
        output_excel_file = '%s/%ipeptides_%iperpool_%ix.xlsx' % (OUTPUT_DIR, num_peptides, num_peptides_per_pool, num_coverage)
        df_design = block_design.to_dataframe()
        writer = pd.ExcelWriter(output_excel_file, engine='openpyxl')
        df_assignment.to_excel(writer, sheet_name='block_assignment', index=False)
        df_design.to_excel(writer, sheet_name='block_design', index=False)
        writer.save()
        print('Finished generating: %i peptides / %i peptides per pool / %ix coverage' %
              (num_peptides, num_peptides_per_pool, num_coverage))

    df_run = pd.DataFrame(data)
    df_run.to_csv(OUTPUT_DIR + '/designs_optimality.csv', index=False)

