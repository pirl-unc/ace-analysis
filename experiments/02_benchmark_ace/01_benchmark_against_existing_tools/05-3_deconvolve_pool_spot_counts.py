import os
import glob
import re
import csv
import logging
import numpy as np
import pandas as pd
import multiprocessing as mp
import time
import shutil
from pathlib import Path
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from acelib.deconvolution import DeconvolutionResult
from acelib.constants import DeconvolutionLabels


PEPTIDES_DIR = "/Users/leework/Documents/Research/projects/project_ace/data/processed/02_benchmark_ace/01_benchmark_against_existing_tools/01_sampled_peptides"
CONFIGURATION_DIR = "/Users/leework/Documents/Research/projects/project_ace/data/processed/02_benchmark_ace/01_benchmark_against_existing_tools/03_generated_configurations/ace-s"
POOL_SPOT_COUNTS_DIR = "/Users/leework/Documents/Research/projects/project_ace/data/processed/02_benchmark_ace/01_benchmark_against_existing_tools/04_simulated_pool_spot_counts/300immunogenic_30nonimmunogenic_1dispersion_fnr0.00/ace-s"
OUTPUT_DIR = "/Users/leework/Documents/Research/projects/project_ace/data/processed/02_benchmark_ace/01_benchmark_against_existing_tools/05_deconvolution_results/300immunogenic_30nonimmunogenic_1dispersion_fnr0.00/strom"
WAIT_TIME = 5
NUM_PROCESSES = 3


def get_logger(name: str) -> logging:
    logging.basicConfig(
        format = '%(asctime)s %(levelname)-8s %(message)s',
        level = logging.INFO,
        datefmt = '%Y-%m-%d %H:%M:%S'
    )
    return logging.getLogger(name)


logger = get_logger(__name__)


def convert_ace_config_to_strom(
        ace_config_file_path,
        write_to_file=False,
        output_file_path='design_file.csv'
):
    """
    Convert the ACE configuration to the form expected by Strom et al 2015

    Parameters:

    ace_config_file_path: Path to the ACE configuration file (.tsv) with the following columns:
        - coverage_id
        - pool_id
        - peptide_id
        - peptide_sequence
        - plate_id
        - well_id

    Returns:

    strom_config_df: The configuration in the form expected by Strom et al 2015 (single column no header)
                     index is the pool_id and the values are comma separated peptide IOs
    """
    df_reduced = pd.read_csv(ace_config_file_path, sep='\t', usecols=['pool_id', 'peptide_id']) # Take only the relevant columns
    pooled_peptides = list(df_reduced.groupby('pool_id')['peptide_id'].apply(lambda x: ', '.join(x))) # Group by well_id and join the peptide_ids
    reformatted_pooled_peptides = [[int(x) for x in re.sub(' ', '', (re.sub('peptide_', '', pooled_peptides[i]))).split(',')] for i in range(len(pooled_peptides))] # Reformat to get rid of the peptide_ and spaces

    # Iterate through the list of lists and convert to a single list
    strom_config = []
    for i in range(len(reformatted_pooled_peptides)):
        strom_config.append(','.join([str(x) for x in reformatted_pooled_peptides[i]]))
    # Convert to a dataframe
    strom_config_df = pd.DataFrame(strom_config, columns=[''])
    if not write_to_file:
        return strom_config_df
    else:
        # Write the list of lists to the CSV file
        with open(output_file_path, mode='w', newline='') as file:
            writer = csv.writer(file, delimiter=',', quoting=csv.QUOTE_NONE)
            writer.writerows(reformatted_pooled_peptides)
    return


def convert_ace_spot_counts_to_strom(
        ace_spot_counts_file_path,
        write_to_file=False,
        output_file_path='pooled_response.csv'
):
    """
    Covnert the ACE spot counts to the form expected by Strom et al 2015

    Parameters:

    ace_spot_counts_file_path: Path to the ACE spot counts file (.tsv) with the following columns:
        - pool_id
        - plate_id
        - well_id
        - spot_count

    """
    df_w_neg = pd.read_csv(ace_spot_counts_file_path, sep='\t')
    # Create a mask for negative wells
    mask = df_w_neg['pool_id'] == 'negative'

    # Get the maximum integer value of the pool)_ids
    df_w_neg.loc[~mask, 'pool_id'] = df_w_neg.loc[~mask, 'pool_id'].astype(int)
    max_pool_id = df_w_neg.loc[~mask, 'pool_id'].max()

    # Fill the negative wells' pool_id with the next integer value
    df_w_neg.loc[mask, 'pool_id'] = np.arange(max_pool_id + 1, max_pool_id + 1 + mask.sum())
    df_w_neg.sort_values(by=['pool_id'], inplace=True)

    if not write_to_file:
        return df_w_neg[['spot_count']].reset_index(drop=True)
    else:
        df_w_neg[['spot_count']].reset_index(drop=True).to_csv(output_file_path, index=False, header=False)
    return


def deconvolve_strom(
        driver,
        design_csv_file: str,
        pool_spot_counts_csv_file: str
):
    # Convert the design file to the form expected by Strom 2015
    try:
        driver.get('https://elispot.shinyapps.io/Shiny/')

        # Upload design file
        design_file = driver.find_element(By.XPATH, '/html/body/div[1]/div/div[1]/div/div[1]/form/div[1]/div[1]/input')
        design_file.send_keys(design_csv_file)

        # Upload pool spot counts file
        elispot_data_file = driver.find_element(By.XPATH, '/html/body/div[1]/div/div[1]/div/div[1]/form/div[1]/div[2]/input')
        elispot_data_file.send_keys(pool_spot_counts_csv_file)
        time.sleep(WAIT_TIME)

        # Click "OK" and wait
        analysis_btn = driver.find_element(By.XPATH, '/html/body/div[1]/div/div[1]/div/div[1]/form/div[3]/button')
        analysis_btn.click()
        time.sleep(WAIT_TIME)

        # Click "Download"
        download_btn = driver.find_element(By.XPATH, '/html/body/div[1]/div/div[1]/div/div[1]/form/a')
        download_btn.click()
        time.sleep(WAIT_TIME)
        return True
    except:
        return False


def process_file(output_dir, intermediate_dir, configurations_dir, pool_spot_counts_csv_files, shared_list):
    chrome_options = Options()
    chrome_options.headless = True
    chrome_options.add_experimental_option("prefs", {
        "download.default_directory": intermediate_dir,
        "download.prompt_for_download": False,
        "download.directory_upgrade": True,
        "safebrowsing.enabled": True
    })
    driver = webdriver.Chrome(options=chrome_options)
    data = {
        'file': [],
        'successful': []
    }
    for pool_spot_counts_csv_file in pool_spot_counts_csv_files:
        pool_spot_counts_csv_file_basename = os.path.basename(pool_spot_counts_csv_file)
        configuration_csv_file = configurations_dir + '/' + pool_spot_counts_csv_file_basename.replace('_pool_spot_counts','')
        output_file = output_dir + '/' + pool_spot_counts_csv_file_basename.replace('_strom_ready.csv', '_strom_results.csv')
        if Path(output_file).is_file():
            logger.info('%s was already downloaded' % output_file)
            data['file'].append(pool_spot_counts_csv_file_basename)
            data['successful'].append(True)
            continue
        is_successful = deconvolve_strom(
            driver=driver,
            design_csv_file=configuration_csv_file,
            pool_spot_counts_csv_file=pool_spot_counts_csv_file
        )
        if is_successful:
            # Rename file
            old_results_csv_file = intermediate_dir + '/Result.csv'
            new_results_csv_file = intermediate_dir + '/' + pool_spot_counts_csv_file_basename.replace('_strom_ready.csv', '_strom_results.csv')
            os.rename(old_results_csv_file, new_results_csv_file)
            shutil.move(new_results_csv_file, output_file)
            logger.info('%s was successfully downloaded' % pool_spot_counts_csv_file_basename)
        else:
            logger.info("%s failed" % pool_spot_counts_csv_file_basename)
        data['file'].append(pool_spot_counts_csv_file_basename)
        data['successful'].append(is_successful)
    df_executions = pd.DataFrame(data)
    shared_list.append(df_executions)


def run_strom_deconvolution(pool_spot_counts_dir, configurations_dir, output_dir):
    pool_spot_counts_csv_files = sorted(glob.glob(pool_spot_counts_dir + '/120*pool_spot_counts_strom_ready.csv'))
    pool_spot_counts_csv_files_list = np.array_split(np.array(pool_spot_counts_csv_files), NUM_PROCESSES)
    pool = mp.Pool(processes=NUM_PROCESSES)
    idx = 1
    manager = mp.Manager()
    L = manager.list()
    for lst in pool_spot_counts_csv_files_list:
        output_dir_ = '%s/intermediate_%i' % (output_dir, idx)
        if not os.path.exists(output_dir_):
            os.makedirs(output_dir_)
        pool.apply_async(process_file, args=[output_dir, output_dir_, configurations_dir, lst, L])
        idx += 1
    pool.close()
    pool.join()
    df_executions = pd.DataFrame()
    for result in L:
        df_executions = pd.concat([df_executions, result], axis=0, ignore_index=True)
    df_executions.to_csv(os.path.join(output_dir, 'execution_logs.tsv'), sep='\t', index=False)


def prepare_strom_configuration_files(configurations_dir, output_dir):
    output_dir_ = output_dir + '/upload_ready'
    if not os.path.exists(output_dir_):
        os.makedirs(output_dir_)
    for design_tsv_file in glob.glob(configurations_dir + '/120*.tsv'):
        design_tsv_file_basename = os.path.basename(design_tsv_file)
        output_file = output_dir_ + '/' + design_tsv_file_basename.replace('.tsv', '_strom_ready.csv')
        convert_ace_config_to_strom(
            ace_config_file_path=design_tsv_file,
            output_file_path=output_file,
            write_to_file=True
        )


def prepare_strom_pool_spot_counts_files(pool_spot_counts_dir, output_dir):
    output_dir_ = output_dir + '/upload_ready'
    for pool_spot_counts_tsv_file in glob.glob(pool_spot_counts_dir + '/120*_pool_spot_counts.tsv'):
        pool_spot_counts_tsv_file_basename = os.path.basename(pool_spot_counts_tsv_file)
        output_file = output_dir_ + '/' + pool_spot_counts_tsv_file_basename.replace('.tsv', '_strom_ready.csv')
        convert_ace_spot_counts_to_strom(
            ace_spot_counts_file_path=pool_spot_counts_tsv_file,
            output_file_path=output_file,
            write_to_file=True
        )


def convert_to_deconvolution_results_worker(peptides_tsv_file, strom_results_csv_file, output_file):
    df_peptides = pd.read_csv(peptides_tsv_file, sep='\t')
    df_strom_results = pd.read_csv(strom_results_csv_file)
    deconvolution_result = DeconvolutionResult()
    for _, row in df_strom_results.iterrows():
        peptide_id = 'peptide_%i' % int(row['Peptide'])
        if peptide_id in df_peptides['peptide_id'].unique():
            peptide_sequence = df_peptides.loc[df_peptides['peptide_id'] == peptide_id, 'epitope'].values.tolist()[0]
            estimated_peptide_spot_count = float(row['  Estimate'])
            if estimated_peptide_spot_count > 0:
                label = DeconvolutionLabels.CANDIDATE_HIT
            else:
                label = DeconvolutionLabels.NOT_A_HIT
            deconvolution_result.add_peptide(
                peptide_id=peptide_id,
                peptide_sequence=peptide_sequence,
                estimated_peptide_spot_count=estimated_peptide_spot_count,
                label=label,
                hit_pool_ids=[]
            )
    deconvolution_result.to_dataframe().to_csv(output_file, sep='\t', index=False)


def convert_to_deconvolution_results(output_dir):
    pool = mp.Pool(processes=NUM_PROCESSES)
    for strom_results_csv_file in glob.glob(output_dir + '/*strom_results.csv'):
        strom_results_csv_file_basename = os.path.basename(strom_results_csv_file)
        peptides_tsv_file = PEPTIDES_DIR + '/' + '_'.join(strom_results_csv_file_basename.split('_')[0:3]) + '.tsv'
        output_file = output_dir + '/' + strom_results_csv_file_basename.replace('.csv', '.tsv')
        pool.apply_async(convert_to_deconvolution_results_worker, args=[peptides_tsv_file,
                                                                        strom_results_csv_file,
                                                                        output_file])
    pool.close()
    pool.join()


if __name__ == "__main__":
    # prepare_strom_configuration_files(
    #     configurations_dir=CONFIGURATION_DIR,
    #     output_dir=OUTPUT_DIR
    # )
    # prepare_strom_pool_spot_counts_files(
    #     pool_spot_counts_dir=POOL_SPOT_COUNTS_DIR,
    #     output_dir=OUTPUT_DIR
    # )
    # run_strom_deconvolution(
    #     pool_spot_counts_dir=OUTPUT_DIR + '/upload_ready',
    #     configurations_dir=OUTPUT_DIR + '/upload_ready',
    #     output_dir=OUTPUT_DIR
    # )
    convert_to_deconvolution_results(output_dir=OUTPUT_DIR)
