"""
The purpose of this python3 script is to evaluate deconvolution results.
"""


import glob
import os
import pandas as pd
import multiprocessing as mp
from acesim.experiment import Experiment
from acelib.block_assignment import BlockAssignment
from acelib.deconvolution import DeconvolutionResult
from acelib.constants import DeconvolutionLabels


NUM_PROCESSES = 32
SAMPLED_PEPTIDES_DIR = "/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/04_real_world_datasets/cameron_2013/01_sampled_peptides"
CONFIGURATIONS_DIR = "/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/04_real_world_datasets/cameron_2013/03_generated_configurations"
DECONVOLUTION_RESULTS_DIR = "/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/04_real_world_datasets/cameron_2013/05_deconvolution_results/300immunogenic_30nonimmunogenic_1dispersion_fnr0.00"
OUTPUT_DIR = "/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/04_real_world_datasets/cameron_2013/06_evaluation_results/300immunogenic_30nonimmunogenic_1dispersion_fnr0.00"


def evaluate_deconvolution(task):
    deconvolution_results_tsv_file = task[0]
    peptides_tsv_file = task[1]
    configuration_tsv_file = task[2]
    configuration_method = task[3]
    deconvolution_method = task[4]
    deconvolution_results_tsv_file_basename = os.path.basename(deconvolution_results_tsv_file)
    df_deconvolution_results = pd.read_csv(deconvolution_results_tsv_file, sep='\t')
    df_deconvolution_results = df_deconvolution_results.loc[df_deconvolution_results['deconvolution_result'] != DeconvolutionLabels.NOT_A_HIT, :]
    deconvolution_result = DeconvolutionResult()
    for _, row in df_deconvolution_results.iterrows():
        peptide_id = row['peptide_id']
        peptide_sequence = row['peptide_sequence']
        estimated_peptide_spot_count = float(row['estimated_peptide_spot_count'])
        hit_pool_ids = [int(pool_id) for pool_id in row['hit_pool_ids'].split(';')]
        deconvolution_label = row['deconvolution_result']
        deconvolution_result.add_peptide(
            peptide_id=peptide_id,
            peptide_sequence=peptide_sequence,
            estimated_peptide_spot_count=estimated_peptide_spot_count,
            hit_pool_ids=hit_pool_ids,
            label=deconvolution_label
        )
    df_peptides = pd.read_csv(peptides_tsv_file, sep='\t')
    df_assignments = pd.read_csv(configuration_tsv_file, sep='\t')
    block_assignments = BlockAssignment.load_from_dataframe(df_assignments=df_assignments)
    rep_id = deconvolution_results_tsv_file_basename.split('_')[2]
    num_peptides = int(deconvolution_results_tsv_file_basename.split('_')[0].replace('peptides', ''))
    num_peptides_per_pool = int(deconvolution_results_tsv_file_basename.split('_')[3].replace('perpool', ''))
    num_coverage = int(deconvolution_results_tsv_file_basename.split('_')[4].replace('x', ''))
    num_immunogenic_peptides = int(deconvolution_results_tsv_file_basename.split('_')[1].replace('immunogenic', ''))
    experiment_id = '%ipeptides_%iimmunogenic_%s_%iperpool_%ix_%s_%s' % (num_peptides,
                                                                         num_immunogenic_peptides,
                                                                         rep_id,
                                                                         num_peptides_per_pool,
                                                                         num_coverage,
                                                                         configuration_method,
                                                                         deconvolution_method)
    df_evaluation_metrics, df_evaluation_results = Experiment.evaluate_deconvolution_results(
        experiment_id=experiment_id,
        df_peptides=df_peptides,
        block_assignment=block_assignments,
        deconvolution_result=deconvolution_result
    )
    df_evaluation_metrics['num_peptides'] = num_peptides
    df_evaluation_metrics['num_peptides_per_pool'] = num_peptides_per_pool
    df_evaluation_metrics['num_coverage'] = num_coverage
    df_evaluation_metrics['num_immunogenic_peptides'] = num_immunogenic_peptides
    df_evaluation_metrics['rep_id'] = rep_id
    df_evaluation_metrics['configuration_method'] = configuration_method
    df_evaluation_metrics['deconvolution_method'] = deconvolution_method
    df_evaluation_results['num_peptides'] = num_peptides
    df_evaluation_results['num_peptides_per_pool'] = num_peptides_per_pool
    df_evaluation_results['num_coverage'] = num_coverage
    df_evaluation_results['num_immunogenic_peptides'] = num_immunogenic_peptides
    df_evaluation_results['rep_id'] = rep_id
    df_evaluation_results['configuration_method'] = configuration_method
    df_evaluation_results['deconvolution_method'] = deconvolution_method
    return (df_evaluation_metrics, df_evaluation_results)


def evaluate(
        deconvolution_results_dir,
        configurations_dir,
        peptides_dir,
        configuration_method,
        deconvolution_method,
        evaluation_metrics_output_file,
        evaluation_results_output_file
):
    tasks = []
    for deconvolution_results_tsv_file in glob.glob(deconvolution_results_dir + '/*pool_spot_counts*.tsv'):
        deconvolution_results_tsv_file_basename = os.path.basename(deconvolution_results_tsv_file)
        peptides_tsv_file = peptides_dir + '/' + '_'.join(deconvolution_results_tsv_file_basename.split('_')[0:3]) + '.tsv'
        configuration_tsv_file = glob.glob(configurations_dir + '/' + '_'.join(deconvolution_results_tsv_file_basename.split('_')[0:5]) + '*.tsv')[0]
        tasks.append((deconvolution_results_tsv_file,
                      peptides_tsv_file,
                      configuration_tsv_file,
                      configuration_method,
                      deconvolution_method))
    pool = mp.Pool(processes=NUM_PROCESSES)
    results = pool.map(evaluate_deconvolution, tasks)
    pool.close()
    pool.join()
    df_evaluation_metrics_all = pd.DataFrame()
    df_evaluation_results_all = pd.DataFrame()
    for result in results:
        df_evaluation_metrics_all = pd.concat([df_evaluation_metrics_all, result[0]], axis=0, ignore_index=True)
        df_evaluation_results_all = pd.concat([df_evaluation_results_all, result[1]], axis=0, ignore_index=True)
    df_evaluation_metrics_all.to_csv(evaluation_metrics_output_file, sep='\t', index=False)
    df_evaluation_results_all.to_csv(evaluation_results_output_file, sep='\t', index=False)


if __name__ == '__main__':
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    # Configuration method  :   ACE-S
    # Deconvolution method  :   ACE-CEM
    evaluate(
        deconvolution_results_dir=DECONVOLUTION_RESULTS_DIR + '/ace-s/ace_cem',
        configurations_dir=CONFIGURATIONS_DIR + '/ace-s',
        peptides_dir=SAMPLED_PEPTIDES_DIR,
        configuration_method='ace-s',
        deconvolution_method='ace-cem',
        evaluation_metrics_output_file=OUTPUT_DIR + '/configuration_ace-s_deconvolution_ace-cem_evaluation_metrics.tsv',
        evaluation_results_output_file=OUTPUT_DIR + '/configuration_ace-s_deconvolution_ace-cem_evaluation_results.tsv'
    )

    # Configuration method  :   ACE
    # Deconvolution method  :   ACE-CEM
    evaluate(
        deconvolution_results_dir=DECONVOLUTION_RESULTS_DIR + '/ace/ace_cem',
        configurations_dir=CONFIGURATIONS_DIR + '/ace',
        peptides_dir=SAMPLED_PEPTIDES_DIR,
        configuration_method='ace',
        deconvolution_method='ace-cem',
        evaluation_metrics_output_file=OUTPUT_DIR + '/configuration_ace_deconvolution_ace-cem_evaluation_metrics.tsv',
        evaluation_results_output_file=OUTPUT_DIR + '/configuration_ace_deconvolution_ace-cem_evaluation_results.tsv'
    )

    # Configuration method  :   Repeated
    # Deconvolution method  :   ACE-EMPIRICAL
    evaluate(
        deconvolution_results_dir=DECONVOLUTION_RESULTS_DIR + '/repeated/ace_empirical',
        configurations_dir=CONFIGURATIONS_DIR + '/repeated',
        peptides_dir=SAMPLED_PEPTIDES_DIR,
        configuration_method='repeated',
        deconvolution_method='ace-empirical',
        evaluation_metrics_output_file=OUTPUT_DIR + '/configuration_repeated_deconvolution_ace-empirical_evaluation_metrics.tsv',
        evaluation_results_output_file=OUTPUT_DIR + '/configuration_repeated_deconvolution_ace-empirical_evaluation_results.tsv'
    )

    # Configuration method  :   Randomized
    # Deconvolution method  :   ACE-CEM
    evaluate(
        deconvolution_results_dir=DECONVOLUTION_RESULTS_DIR + '/randomized/ace_cem',
        configurations_dir=CONFIGURATIONS_DIR + '/randomized',
        peptides_dir=SAMPLED_PEPTIDES_DIR,
        configuration_method='randomized',
        deconvolution_method='ace-cem',
        evaluation_metrics_output_file=OUTPUT_DIR + '/configuration_randomized_deconvolution_ace-cem_evaluation_metrics.tsv',
        evaluation_results_output_file=OUTPUT_DIR + '/configuration_randomized_deconvolution_ace-cem_evaluation_results.tsv'
    )
