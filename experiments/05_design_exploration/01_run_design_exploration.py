import pandas as pd
import numpy as np
import multiprocessing as mp
import os
from acesim.experiment import Experiment
from acesim.configurator_ace import AceConfigurator
from acesim.deconvolver_ace import AceDeconvolver
from acesim.utilities import *
from acelib.deconvolution import DeconvolutionResult
from acelib.constants import DeconvolutionLabels


DESIGN_CONFIGURATIONS_TSV_FILE = "/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/raw/design_exploration/design_configurations.tsv"
REFERENCE_CSV_FILE = "/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/raw/references/large_dummy_sequence_dataset.csv"
OUTPUT_DIR = "/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/05_design_exploration"
NUM_PROCESSES = 96
PEPTIDE_LENGTH = 9
SAMPLING_METHOD = 'random'
NUM_REPLICATES = 100
MU_IMMUNOGENIC = 300.0
MU_NONIMMUNOGENIC = 30.0
DISPERSION_FACTOR = 1.0
FALSE_NEGATIVE_RATE = 0.0
MODE = 'golfy'
SEQUENCE_SIMILARITY_FUNCTION = 'euclidean'
SEQUENCE_SIMILARITY_THRESHOLD = 0.7
GOLFY_STRATEGY = 'greedy'
GOLFY_MAX_ITERS = 2000
PLATE_SIZE = 96
NUM_NEGATIVE_CONTROL_WELLS = 3
NEGATIVE_CONTROL_MULTIPLIER = 3
DECONVOLUTION_METHOD = 'cem'
CONFIGURATION_METHOD = 'ace-s'


def run_experiment(task):
    df_ref_peptides = task[0]
    num_peptides = task[1]
    num_peptides_per_pool = task[2]
    num_coverage = task[3]
    num_peptides_immunogenic = task[4]
    rep_id = task[5]

    # Step 1. Sample peptides
    df_peptides = Experiment.sample_peptides(
        df_ref_peptides=df_ref_peptides,
        num_peptides=num_peptides,
        num_peptides_immunogenic=num_peptides_immunogenic,
        peptide_length=PEPTIDE_LENGTH,
        sampling_method=SAMPLING_METHOD
    )

    # Step 2. Sample peptides spot counts
    df_peptides_spot_counts = Experiment.simulate_peptide_spot_counts(
        df_peptides=df_peptides,
        num_coverage=num_coverage,
        random_effects=True,
        mu_immunogenic=MU_IMMUNOGENIC,
        mu_nonimmunogenic=MU_NONIMMUNOGENIC,
        dispersion_factor=DISPERSION_FACTOR
    )

    # Step 3. Generate configuration
    ace_configurator = AceConfigurator(name='ace-s')
    random_seed = Experiment.generate_random_seed()
    peptides = convert_peptides_dataframe_to_peptides(df_peptides=df_peptides)
    block_assignment, _ = ace_configurator.generate_assignment(
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

    # Step 4. Simulate pool spot counts
    df_assignments = block_assignment.to_dataframe()
    df_readout_ = Experiment.aggregate_pool_spot_counts(
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
    for _, row in df_readout_.iterrows():
        df_matched = df_assignments.loc[df_assignments['pool_id'] == row['pool_id'], :]
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
    df_readout = pd.DataFrame(data)

    # Step 5. Deconvolve
    min_pool_spot_count = np.mean(df_readout.loc[df_readout['pool_id'] == 'negative','spot_count'].values.tolist()) * NEGATIVE_CONTROL_MULTIPLIER
    df_readout = df_readout.loc[df_readout['pool_id'] != 'negative',:]
    df_readout['pool_id'] = df_readout['pool_id'].astype('int')
    ace_deconvolver = AceDeconvolver(name='ace')
    deconvolution_result = ace_deconvolver.deconvolve(
        df_readout=df_readout,
        block_assignment=block_assignment,
        method=DECONVOLUTION_METHOD,
        min_coverage=num_coverage,
        min_pool_spot_count=min_pool_spot_count
    )
    df_deconvolution_results = deconvolution_result.to_dataframe()

    # Step 6. Evaluate
    df_deconvolution_results = df_deconvolution_results.loc[df_deconvolution_results['deconvolution_result'] != DeconvolutionLabels.NOT_A_HIT, :]
    deconvolution_result_ = DeconvolutionResult()
    for _, row in df_deconvolution_results.iterrows():
        peptide_id = row['peptide_id']
        peptide_sequence = row['peptide_sequence']
        estimated_peptide_spot_count = float(row['estimated_peptide_spot_count'])
        hit_pool_ids = [int(pool_id) for pool_id in row['hit_pool_ids'].split(';')]
        deconvolution_label = row['deconvolution_result']
        deconvolution_result_.add_peptide(
            peptide_id=peptide_id,
            peptide_sequence=peptide_sequence,
            estimated_peptide_spot_count=estimated_peptide_spot_count,
            hit_pool_ids=hit_pool_ids,
            label=deconvolution_label
        )
    experiment_id = '%ipeptides_%iimmunogenic_%s_%iperpool_%ix_%s_%s' % (
        num_peptides,
        num_peptides_immunogenic,
        rep_id,
        num_peptides_per_pool,
        num_coverage,
        CONFIGURATION_METHOD,
        DECONVOLUTION_METHOD
    )
    df_evaluation_metrics, df_evaluation_results = Experiment.evaluate_deconvolution_results(
        experiment_id=experiment_id,
        df_peptides=df_peptides,
        block_assignment=block_assignment,
        deconvolution_result=deconvolution_result_
    )
    df_evaluation_metrics['num_peptides'] = num_peptides
    df_evaluation_metrics['num_peptides_per_pool'] = num_peptides_per_pool
    df_evaluation_metrics['num_coverage'] = num_coverage
    df_evaluation_metrics['num_immunogenic_peptides'] = num_peptides_immunogenic
    df_evaluation_metrics['rep_id'] = rep_id
    df_evaluation_metrics['configuration_method'] = CONFIGURATION_METHOD
    df_evaluation_metrics['deconvolution_method'] = DECONVOLUTION_METHOD
    return df_evaluation_metrics


if __name__ == "__main__":
    print("Started")
    df_designs = pd.read_csv(DESIGN_CONFIGURATIONS_TSV_FILE, sep='\t')
    df_ref_peptides = pd.read_csv(REFERENCE_CSV_FILE)
    tasks = []
    for _, row in df_designs.iterrows():
        num_peptides = row['num_peptides']
        num_peptides_per_pool = row['num_peptides_per_pool']
        num_coverage = row['num_coverage']
        num_peptides_immunogenic = row['num_immunogenic_peptides']
        for i in range(0, NUM_REPLICATES):
            tasks.append((
                df_ref_peptides,
                num_peptides,
                num_peptides_per_pool,
                num_coverage,
                num_peptides_immunogenic,
                'rep%i' % (i + 1)
            ))
    pool = mp.Pool(processes=NUM_PROCESSES)
    results = pool.map(run_experiment, tasks)
    pool.close()
    pool.join()
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    df_evaluation_metrics_all = pd.DataFrame()
    for result in results:
        df_evaluation_metrics_all = pd.concat([df_evaluation_metrics_all, result], axis=0, ignore_index=True)
    df_evaluation_metrics_all.to_csv(OUTPUT_DIR + '/evaluations.tsv', sep='\t', index=False)
    print("Finished successfully")
