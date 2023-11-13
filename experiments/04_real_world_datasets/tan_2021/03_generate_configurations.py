"""
The purpose of this python3 script is to generate configurations.
"""


import glob
import pandas as pd
import os
import multiprocessing as mp
from acesim.configurator_ace import AceConfigurator
from acesim.configurator_randomized import RandomizedDesignConfigurator
from acesim.configurator_repeated import RepeatedDesignConfigurator
from acesim.utilities import *
from acesim.experiment import Experiment


NUM_PROCESSES = 8
SAMPLED_PEPTIDES_DIR = "/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/04_real_world_datasets/tan_2021/01_sampled_peptides"
OUTPUT_DIR = "/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/04_real_world_datasets/tan_2021/03_generated_configurations"
TRAINED_MODEL_FILE = "/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/raw/trained_models/trained_model_w_data_augmentation_b3000.pt"
MODE = 'golfy'
SEQUENCE_SIMILARITY_FUNCTION = 'euclidean'
SEQUENCE_SIMILARITY_THRESHOLD = 0.7
GOLFY_STRATEGY = 'greedy'
GOLFY_MAX_ITERS = 2000
PLATE_SIZE = 96


def generate_repeated_design(task):
    peptide_tsv_file = task[0]
    num_peptides = task[1]
    num_peptides_per_pool = task[2]
    num_coverage = task[3]
    file_basename = os.path.basename(peptide_tsv_file)
    df_peptides = pd.read_csv(peptide_tsv_file, sep='\t')
    peptides = convert_peptides_dataframe_to_peptides(df_peptides=df_peptides)
    repeated_configurator = RepeatedDesignConfigurator(name='repeated')
    block_assignment, _ = repeated_configurator.generate_assignment(
        peptides=peptides,
        num_peptides_per_pool=num_peptides_per_pool,
        num_coverage=num_coverage,
        plate_size=PLATE_SIZE
    )
    output_file = '%s_%iperpool_%ix_repeated_design.tsv' % (file_basename.replace('.tsv', ''), num_peptides_per_pool, num_coverage)
    block_assignment.to_dataframe().to_csv(OUTPUT_DIR + '/repeated/' + output_file,
                                           sep='\t', index=False)


def generate_randomized_design(task):
    peptide_tsv_file = task[0]
    num_peptides = task[1]
    num_peptides_per_pool = task[2]
    num_coverage = task[3]
    file_basename = os.path.basename(peptide_tsv_file)
    df_peptides = pd.read_csv(peptide_tsv_file, sep='\t')
    peptides = convert_peptides_dataframe_to_peptides(df_peptides=df_peptides)
    randomized_configurator = RandomizedDesignConfigurator(name='randomized')
    random_seed = Experiment.generate_random_seed()
    block_assignment, _ = randomized_configurator.generate_assignment(
        peptides=peptides,
        num_peptides_per_pool=num_peptides_per_pool,
        num_coverage=num_coverage,
        random_seed=random_seed,
        plate_size=PLATE_SIZE
    )
    output_file = '%s_%iperpool_%ix_randomized_design_%i.tsv' % (file_basename.replace('.tsv', ''), num_peptides_per_pool, num_coverage, random_seed)
    block_assignment.to_dataframe().to_csv(OUTPUT_DIR + '/randomized/' + output_file,
                                           sep='\t', index=False)


def generate_ace_design_without_sequence_similarity(task):
    peptide_tsv_file = task[0]
    num_peptides = task[1]
    num_peptides_per_pool = task[2]
    num_coverage = task[3]
    file_basename = os.path.basename(peptide_tsv_file)
    df_peptides = pd.read_csv(peptide_tsv_file, sep='\t')
    peptides = convert_peptides_dataframe_to_peptides(df_peptides=df_peptides)
    ace_configurator = AceConfigurator(name='ace')
    random_seed = Experiment.generate_random_seed()
    block_assignment, preferred_peptide_pairs = ace_configurator.generate_assignment(
        peptides=peptides,
        num_peptides_per_pool=num_peptides_per_pool,
        num_coverage=num_coverage,
        trained_model_file='',
        cluster_peptides=False,
        mode=MODE,
        sequence_similarity_function=SEQUENCE_SIMILARITY_FUNCTION,
        sequence_similarity_threshold=SEQUENCE_SIMILARITY_THRESHOLD,
        golfy_random_seed=random_seed,
        golfy_strategy=GOLFY_STRATEGY,
        golfy_max_iters=GOLFY_MAX_ITERS,
        golfy_allow_extra_pools=False,
        plate_size=PLATE_SIZE
    )

    # Block assignment
    output_file = '%s_%iperpool_%ix_ace_design_%i.tsv' % (file_basename.replace('.tsv', ''), num_peptides_per_pool, num_coverage, random_seed)
    block_assignment.to_dataframe().to_csv(OUTPUT_DIR + '/ace/' + output_file,
                                           sep='\t', index=False)


def generate_ace_design(task):
    peptide_tsv_file = task[0]
    num_peptides = task[1]
    num_peptides_per_pool = task[2]
    num_coverage = task[3]
    file_basename = os.path.basename(peptide_tsv_file)
    df_peptides = pd.read_csv(peptide_tsv_file, sep='\t')
    peptides = convert_peptides_dataframe_to_peptides(df_peptides=df_peptides)
    ace_configurator = AceConfigurator(name='ace-s')
    random_seed = Experiment.generate_random_seed()
    block_assignment, preferred_peptide_pairs = ace_configurator.generate_assignment(
        peptides=peptides,
        num_peptides_per_pool=num_peptides_per_pool,
        num_coverage=num_coverage,
        trained_model_file=TRAINED_MODEL_FILE,
        cluster_peptides=True,
        mode=MODE,
        sequence_similarity_function=SEQUENCE_SIMILARITY_FUNCTION,
        sequence_similarity_threshold=SEQUENCE_SIMILARITY_THRESHOLD,
        golfy_random_seed=random_seed,
        golfy_strategy=GOLFY_STRATEGY,
        golfy_max_iters=GOLFY_MAX_ITERS,
        golfy_allow_extra_pools=False,
        plate_size=PLATE_SIZE
    )

    # Block assignment
    output_file = '%s_%iperpool_%ix_ace-s_design_%i.tsv' % (file_basename.replace('.tsv', ''), num_peptides_per_pool, num_coverage, random_seed)
    block_assignment.to_dataframe().to_csv(OUTPUT_DIR + '/ace-s/' + output_file,
                                           sep='\t', index=False)

    # Preferred peptide pairs
    data = {
        'peptide_id_1': [],
        'peptide_id_2': [],
        'peptide_sequence_1': [],
        'peptide_sequence_2': [],
        'similarity_score': []
    }
    for peptide_id_1, peptide_id_2, similarity_score in preferred_peptide_pairs:
        peptide_sequence_1 = df_peptides.loc[df_peptides['peptide_id'] == peptide_id_1, 'epitope'].values.tolist()[0]
        peptide_sequence_2 = df_peptides.loc[df_peptides['peptide_id'] == peptide_id_2, 'epitope'].values.tolist()[0]
        data['peptide_id_1'].append(peptide_id_1)
        data['peptide_id_2'].append(peptide_id_2)
        data['peptide_sequence_1'].append(peptide_sequence_1)
        data['peptide_sequence_2'].append(peptide_sequence_2)
        data['similarity_score'].append(similarity_score)
    df_preferred_pairs = pd.DataFrame(data)
    output_file = '%s_%iperpool_%ix_ace-s_design_%i_preferred_peptide_pairs.tsv' % (file_basename.replace('.tsv', ''), num_peptides_per_pool, num_coverage, random_seed)
    df_preferred_pairs.to_csv(OUTPUT_DIR + '/ace-s/preferred_peptide_pairs/' + output_file,
                              sep='\t', index=False)


if __name__ == "__main__":
    if not os.path.exists(OUTPUT_DIR + '/ace-s/'):
        os.makedirs(OUTPUT_DIR + '/ace-s/')
    if not os.path.exists(OUTPUT_DIR + '/ace-s/preferred_peptide_pairs/'):
        os.makedirs(OUTPUT_DIR + '/ace-s/preferred_peptide_pairs/')
    if not os.path.exists(OUTPUT_DIR + '/ace/'):
        os.makedirs(OUTPUT_DIR + '/ace/')
    if not os.path.exists(OUTPUT_DIR + '/repeated/'):
        os.makedirs(OUTPUT_DIR + '/repeated/')
    if not os.path.exists(OUTPUT_DIR + '/randomized/'):
        os.makedirs(OUTPUT_DIR + '/randomized/')

    # Step 1. Create tasks
    tasks = []
    for peptide_tsv_file in glob.glob('%s/*immunogenic_rep*.tsv' % (SAMPLED_PEPTIDES_DIR)):
        tasks.append((peptide_tsv_file, 177, 10, 3))

    # Step 2. Multiprocess tasks
    pool = mp.Pool(processes=NUM_PROCESSES)

    # Generate repeated design configurations
    pool.map(generate_repeated_design, tasks)

    # Generate randomized design configurations
    pool.map(generate_randomized_design, tasks)

    # Generate ACE-S design configurations
    pool.map(generate_ace_design, tasks)

    # Generate ACE design configurations
    pool.map(generate_ace_design_without_sequence_similarity, tasks)

    pool.close()
    pool.join()
