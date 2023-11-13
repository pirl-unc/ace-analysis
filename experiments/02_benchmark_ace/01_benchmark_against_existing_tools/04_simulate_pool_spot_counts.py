"""
The purpose of this python3 script is to simulate pool spot counts.
"""


import glob
import openpyxl
import os
import pandas as pd
import numpy as np
import multiprocessing as mp
from acesim.experiment import Experiment


NUM_PROCESSES = 64
DESIGN_CONFIGURATIONS_TSV_FILE = "/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/raw/benchmark_against_existing_tools/design_configurations.tsv"
SAMPLED_PEPTIDES_SPOT_COUNTS_DIR = "/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/02_benchmark_ace/01_benchmark_against_existing_tools/02_sampled_peptides_spot_counts/300immunogenic_30nonimmunogenic_1dispersion"
STRANDBERG_TEMPLATES_DIR = "/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/raw/references"
CONFIGURATIONS_DIR = "/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/02_benchmark_ace/01_benchmark_against_existing_tools/03_generated_configurations"
OUTPUT_DIR = "/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/02_benchmark_ace/01_benchmark_against_existing_tools/04_simulated_pool_spot_counts/300immunogenic_30nonimmunogenic_1dispersion_fnr0.00"
MU_IMMUNOGENIC = 300.0
MU_NONIMMUNOGENIC = 30.0
DISPERSION_FACTOR = 1.0
SATURATED_SPOT_COUNT = 600
FALSE_NEGATIVE_RATE = 0.00
NUM_NEGATIVE_CONTROL_WELLS = 3


def simulate_pool_spot_counts_worker(task):
    design_tsv_file = task[0]
    peptides_spot_counts_tsv_file = task[1]
    output_file = task[2]
    num_peptides_per_pool = task[3]
    df_assignments = pd.read_csv(design_tsv_file, sep='\t')
    df_peptides_spot_counts = pd.read_csv(peptides_spot_counts_tsv_file, sep='\t')
    df_readout = Experiment.aggregate_pool_spot_counts(
        df_assignments=df_assignments,
        df_peptides_spot_counts=df_peptides_spot_counts,
        mu_nonimmunogenic=MU_NONIMMUNOGENIC,
        dispersion_factor=DISPERSION_FACTOR,
        false_negative_rate=FALSE_NEGATIVE_RATE,
        num_peptides_per_pool=num_peptides_per_pool
    )
    data = {
        'pool_id': [],
        'plate_id': [],
        'well_id': [],
        'spot_count': []
    }
    for _, row in df_readout.iterrows():
        df_matched = df_assignments.loc[df_assignments['pool_id'] == row['pool_id'],:]
        data['pool_id'].append(int(row['pool_id']))
        data['spot_count'].append(round(row['spot_count']))
        data['plate_id'].append(df_matched['plate_id'].values.tolist()[0])
        data['well_id'].append(df_matched['well_id'].values.tolist()[0])

    # Simulate negative control well spot counts
    for _ in range(0, NUM_NEGATIVE_CONTROL_WELLS):
        negative_pool_spot_count = Experiment.sample_spot_counts(
            mean=MU_NONIMMUNOGENIC,
            dispersion_factor=DISPERSION_FACTOR,
            num_coverage=num_peptides_per_pool
        )
        data['pool_id'].append('negative')
        data['plate_id'].append(1)
        data['well_id'].append('')
        data['spot_count'].append(round(np.sum(np.array(negative_pool_spot_count) / float(num_peptides_per_pool))))
    df_readout_ = pd.DataFrame(data)
    df_readout_.to_csv(output_file, sep='\t', index=False)


def simulate_pool_spot_counts():
    if not os.path.exists(OUTPUT_DIR + '/ace-s/'):
        os.makedirs(OUTPUT_DIR + '/ace-s/')
    if not os.path.exists(OUTPUT_DIR + '/ace/'):
        os.makedirs(OUTPUT_DIR + '/ace/')
    if not os.path.exists(OUTPUT_DIR + '/repeated/'):
        os.makedirs(OUTPUT_DIR + '/repeated/')
    if not os.path.exists(OUTPUT_DIR + '/randomized/'):
        os.makedirs(OUTPUT_DIR + '/randomized/')
    if not os.path.exists(OUTPUT_DIR + '/strandberg/'):
        os.makedirs(OUTPUT_DIR + '/strandberg/')

    # Step 1. Create tasks
    df_designs = pd.read_csv(DESIGN_CONFIGURATIONS_TSV_FILE, sep='\t')
    tasks = []
    for _, row in df_designs.iterrows():
        num_peptides = row['num_peptides']
        num_peptides_per_pool = row['num_peptides_per_pool']
        num_coverage = row['num_coverage']
        num_immunogenic_peptides = row['num_immunogenic_peptides']

        # Repeated design
        for design_tsv_file in glob.glob('%s/repeated/%ipeptides_%iimmunogenic_rep*_%iperpool_%ix_repeated_design.tsv' % (CONFIGURATIONS_DIR, num_peptides, num_immunogenic_peptides, num_peptides_per_pool, num_coverage)):
            design_tsv_file_basename = os.path.basename(design_tsv_file)
            design_tsv_file_basename_elements = design_tsv_file_basename.split('_')
            replicate = design_tsv_file_basename_elements[2]
            peptides_spots_count_file = '%s/%ipeptides_%iimmunogenic_%s_%ix_spot_counts.tsv' % (SAMPLED_PEPTIDES_SPOT_COUNTS_DIR, num_peptides, num_immunogenic_peptides, replicate, num_coverage)
            output_file = '%s/repeated/%s' % (OUTPUT_DIR, design_tsv_file_basename.replace('.tsv', '_pool_spot_counts.tsv'))
            tasks.append((design_tsv_file, peptides_spots_count_file, output_file, num_peptides_per_pool))

        # Randomized design
        for design_tsv_file in glob.glob('%s/randomized/%ipeptides_%iimmunogenic_rep*_%iperpool_%ix_randomized_design_*.tsv' % (CONFIGURATIONS_DIR, num_peptides, num_immunogenic_peptides, num_peptides_per_pool, num_coverage)):
            design_tsv_file_basename = os.path.basename(design_tsv_file)
            design_tsv_file_basename_elements = design_tsv_file_basename.split('_')
            replicate = design_tsv_file_basename_elements[2]
            peptides_spots_count_file = '%s/%ipeptides_%iimmunogenic_%s_%ix_spot_counts.tsv' % (SAMPLED_PEPTIDES_SPOT_COUNTS_DIR, num_peptides, num_immunogenic_peptides, replicate, num_coverage)
            output_file = '%s/randomized/%s' % (OUTPUT_DIR, design_tsv_file_basename.replace('.tsv', '_pool_spot_counts.tsv'))
            tasks.append((design_tsv_file, peptides_spots_count_file, output_file, num_peptides_per_pool))

        # ACE-S design
        for design_tsv_file in glob.glob('%s/ace-s/%ipeptides_%iimmunogenic_rep*_%iperpool_%ix_ace-s_design_*.tsv' % (CONFIGURATIONS_DIR, num_peptides, num_immunogenic_peptides, num_peptides_per_pool, num_coverage)):
            design_tsv_file_basename = os.path.basename(design_tsv_file)
            design_tsv_file_basename_elements = design_tsv_file_basename.split('_')
            replicate = design_tsv_file_basename_elements[2]
            peptides_spots_count_file = '%s/%ipeptides_%iimmunogenic_%s_%ix_spot_counts.tsv' % (SAMPLED_PEPTIDES_SPOT_COUNTS_DIR, num_peptides, num_immunogenic_peptides, replicate, num_coverage)
            output_file = '%s/ace-s/%s' % (OUTPUT_DIR, design_tsv_file_basename.replace('.tsv', '_pool_spot_counts.tsv'))
            tasks.append((design_tsv_file, peptides_spots_count_file, output_file, num_peptides_per_pool))

        # ACE design
        for design_tsv_file in glob.glob('%s/ace/%ipeptides_%iimmunogenic_rep*_%iperpool_%ix_ace_design_*.tsv' % (CONFIGURATIONS_DIR, num_peptides, num_immunogenic_peptides, num_peptides_per_pool, num_coverage)):
            design_tsv_file_basename = os.path.basename(design_tsv_file)
            design_tsv_file_basename_elements = design_tsv_file_basename.split('_')
            replicate = design_tsv_file_basename_elements[2]
            peptides_spots_count_file = '%s/%ipeptides_%iimmunogenic_%s_%ix_spot_counts.tsv' % (SAMPLED_PEPTIDES_SPOT_COUNTS_DIR, num_peptides, num_immunogenic_peptides, replicate, num_coverage)
            output_file = '%s/ace/%s' % (OUTPUT_DIR, design_tsv_file_basename.replace('.tsv', '_pool_spot_counts.tsv'))
            tasks.append((design_tsv_file, peptides_spots_count_file, output_file, num_peptides_per_pool))

        # Strandberg design
        if num_peptides == 120:
            num_peptides_per_pool = 11
        if num_peptides == 800:
            num_peptides_per_pool = 28
        for design_tsv_file in glob.glob('%s/strandberg/%ipeptides_%iimmunogenic_rep*_%iperpool_%ix_strandberg_design.tsv' % (CONFIGURATIONS_DIR, num_peptides, num_immunogenic_peptides, num_peptides_per_pool, num_coverage)):
            design_tsv_file_basename = os.path.basename(design_tsv_file)
            design_tsv_file_basename_elements = design_tsv_file_basename.split('_')
            replicate = design_tsv_file_basename_elements[2]
            peptides_spots_count_file = '%s/%ipeptides_%iimmunogenic_%s_%ix_spot_counts.tsv' % (SAMPLED_PEPTIDES_SPOT_COUNTS_DIR, num_peptides, num_immunogenic_peptides, replicate, num_coverage)
            output_file = '%s/strandberg/%s' % (OUTPUT_DIR, design_tsv_file_basename.replace('.tsv', '_pool_spot_counts.tsv'))
            tasks.append((design_tsv_file, peptides_spots_count_file, output_file, num_peptides_per_pool))

    # Step 2. Multiprocess tasks
    pool = mp.Pool(processes=NUM_PROCESSES)
    pool.map(simulate_pool_spot_counts_worker, tasks)
    pool.close()
    pool.join()


def convert_pool_spot_counts_to_strandberg_tool_format_worker(pool_spot_counts_tsv_file):
    pool_spot_counts_tsv_file_basename = os.path.basename(pool_spot_counts_tsv_file)
    num_peptides = int(pool_spot_counts_tsv_file_basename.split('_')[0].replace('peptides', ''))
    if num_peptides == 120:
        template_csv_file = STRANDBERG_TEMPLATES_DIR + '/120peptides_11perpool_3x_strandberg_template.csv'
        num_peptides_per_pool = 11
    if num_peptides == 800:
        template_csv_file = STRANDBERG_TEMPLATES_DIR + '/800peptides_28perpool_3x_strandberg_template.csv'
        num_peptides_per_pool = 28
    df_template = pd.read_csv(template_csv_file)
    df_pool_spot_counts = pd.read_csv(pool_spot_counts_tsv_file, sep='\t')

    # Remove negative wells
    df_pool_spot_counts = df_pool_spot_counts.loc[df_pool_spot_counts['pool_id'] != 'negative',:]
    data = [
        ['Plate:, Date,  User: default, Path:', '', '', '', '', '', '', '', '', '', '', '', ''],
        ['Plate Data - Spots Number - Default count settings', '', '', '', '', '', '', '', '', '', '', '', ''],
        ['', 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
        ['A'],
        ['B'],
        ['C'],
        ['D'],
        ['E'],
        ['F'],
        ['G'],
        ['H']
    ]
    row_idx = 3
    for i in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']:
        for j in range(1, 13):
            well_id = '%s%i' % (i, j)
            df_matched = df_pool_spot_counts.loc[df_pool_spot_counts['well_id'] == well_id, :]
            if len(df_matched) == 0:
                # Add positive and negative control wells
                pool_type = df_template.loc[df_template['Position'] == well_id, 'Antigen'].values.tolist()[0]
                if pool_type == 'Neg':
                    pool_spot_count = Experiment.sample_spot_counts(
                        mean=MU_NONIMMUNOGENIC,
                        dispersion_factor=DISPERSION_FACTOR,
                        num_coverage=num_peptides_per_pool
                    )
                    pool_spot_count = round(np.sum(np.array(pool_spot_count) / float(num_peptides_per_pool)))
                if pool_type == 'Pos':
                    pool_spot_count = Experiment.sample_spot_counts(
                        mean=MU_IMMUNOGENIC,
                        dispersion_factor=DISPERSION_FACTOR,
                        num_coverage=num_peptides_per_pool
                    )
                    pool_spot_count = round(np.sum(pool_spot_count))
                    pool_spot_count = min(round(pool_spot_count), SATURATED_SPOT_COUNT)
            else:
                pool_spot_count = df_matched['spot_count'].values.tolist()[0]
            data[row_idx].append(pool_spot_count)
        row_idx += 1

    # Add the data to the spreadsheet
    workbook = openpyxl.Workbook()
    sheet = workbook.active
    sheet.title = 'PlateData'
    for row_idx, row_data in enumerate(data, start=1):
        for col_idx, value in enumerate(row_data, start=1):
            cell = sheet.cell(row=row_idx, column=col_idx)
            cell.value = value

    # Save the spreadsheet
    output_file = OUTPUT_DIR + '/strandberg/upload_ready/' + pool_spot_counts_tsv_file_basename.replace('.tsv','.xlsx')
    workbook.save(output_file)


def convert_pool_spot_counts_to_strandberg_tool_format():
    pool = mp.Pool(processes=NUM_PROCESSES)
    if not os.path.exists(OUTPUT_DIR + '/strandberg/upload_ready/'):
        os.makedirs(OUTPUT_DIR + '/strandberg/upload_ready/')
    pool.map(convert_pool_spot_counts_to_strandberg_tool_format_worker, glob.glob(OUTPUT_DIR + '/strandberg/*.tsv'))
    pool.close()
    pool.join()


if __name__ == "__main__":
    print("Started")
    simulate_pool_spot_counts()
    convert_pool_spot_counts_to_strandberg_tool_format()
    print("Finished successfully")
