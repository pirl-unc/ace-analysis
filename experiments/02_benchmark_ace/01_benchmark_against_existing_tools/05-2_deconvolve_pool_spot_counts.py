"""
The purpose of this python3 script is to deconvolve pool spot counts using Strandberg's heuristic.
"""


import csv
import sys
import glob
import pandas as pd
import multiprocessing as mp
import logging
import os
import openpyxl
import time
from pathlib import Path
from typing import Literal
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.ui import Select
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support.wait import WebDriverWait
from acelib.deconvolution import DeconvolutionResult
from acelib.constants import DeconvolutionLabels


NUM_PROCESSES = 3
PEPTIDES_DIR = "/Users/leework/Documents/Research/projects/project_ace/data/processed/02_benchmark_ace/01_benchmark_against_existing_tools/01_sampled_peptides"
CONFIGURATIONS_DIR = "/Users/leework/Documents/Research/projects/project_ace/data/processed/02_benchmark_ace/01_benchmark_against_existing_tools/03_generated_configurations/strandberg"
POOL_SPOT_COUNTS_DIR = "/Users/leework/Documents/Research/projects/project_ace/data/processed/02_benchmark_ace/01_benchmark_against_existing_tools/04_simulated_pool_spot_counts"
STRANDBERG_TEMPLATES_DIR = "/Users/leework/Documents/Research/projects/project_ace/data/raw/references"
OUTPUT_DIR = "/Users/leework/Documents/Research/projects/project_ace/data/processed/02_benchmark_ace/01_benchmark_against_existing_tools/05_deconvolution_results"
WAIT_TIME = 5


def get_logger(name: str) -> logging:
    logging.basicConfig(
        format = '%(asctime)s %(levelname)-8s %(message)s',
        level = logging.INFO,
        datefmt = '%Y-%m-%d %H:%M:%S'
    )
    return logging.getLogger(name)


logger = get_logger(__name__)


def deconvolve_strandberg_heuristic(
        driver,
        pool_spot_counts_xlsx_file: str,
        design_csv_file: str,
        positivity_criterion: Literal['background_subtracted', 'empirical', 'bayesian']
):
    try:
        driver.get('https://rickardstrandberg.shinyapps.io/elispotapp/')
        analysis_btn = driver.find_element(By.XPATH, '/html/body/nav/div/ul/li[4]/a')
        analysis_btn.click()

        # Upload pool spot counts file
        elispot_data_file = driver.find_element(By.XPATH, '/html/body/div[1]/div/div[4]/div/div[1]/form/div[1]/div[1]/label/span/input')
        elispot_data_file.send_keys(pool_spot_counts_xlsx_file)

        # Upload design file
        design_file = driver.find_element(By.XPATH, '/html/body/div[1]/div/div[4]/div/div[1]/form/div[2]/div[1]/label/span/input')
        design_file.send_keys(design_csv_file)

        # Set the positivity criterion
        dropdown = driver.find_element(By.XPATH, '/html/body/div[1]/div/div[4]/div/div[1]/form/div[4]/div/div')
        dropdown.click()
        if positivity_criterion == 'background_subtracted':
            option = driver.find_element(By.XPATH, '/html/body/div[1]/div/div[4]/div/div[1]/form/div[4]/div/div/div[2]/div/div[1]')
            option.click()
        elif positivity_criterion == 'empirical':
            option = driver.find_element(By.XPATH, '/html/body/div[1]/div/div[4]/div/div[1]/form/div[4]/div/div/div[2]/div/div[2]')
            option.click()
        elif positivity_criterion == 'bayesian':
            option = driver.find_element(By.XPATH, '/html/body/div[1]/div/div[4]/div/div[1]/form/div[4]/div/div/div[2]/div/div[3]')
            option.click()
        else:
            raise Exception('Unknown positivity criterion: %s' % positivity_criterion)
        time.sleep(WAIT_TIME)

        # Click "Analyze"
        analyse_btn = driver.find_element(By.XPATH, '/html/body/div[1]/div/div[4]/div/div[1]/form/button')
        analyse_btn.click()
        time.sleep(WAIT_TIME)

        # Click "Download"
        download_btn = driver.find_element(By.XPATH, '/html/body/div[1]/div/div[4]/div/div[2]/div[1]/a')
        download_btn.click()
        time.sleep(WAIT_TIME)

        return True
    except:
        return False


def run_strandberg_heuristic(
        pool_spot_counts_dir: str,
        output_dir: str,
        positivity_criterion: str,
        output_file_suffix: str
):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    chrome_options = Options()
    chrome_options.headless = True
    chrome_options.add_experimental_option("prefs", {
        "download.default_directory": output_dir,
        "download.prompt_for_download": False,
        "download.directory_upgrade": True,
        "safebrowsing.enabled": True
    })
    data = {
        'file': [],
        'successful': []
    }
    driver = webdriver.Chrome(options=chrome_options)
    for pool_spot_counts_xlsx_file in sorted(glob.glob(pool_spot_counts_dir + '/*.xlsx')):
        pool_spot_counts_xlsx_file_basename = os.path.basename(pool_spot_counts_xlsx_file)
        new_elispot_design_xlsx_file = output_dir + '/' + pool_spot_counts_xlsx_file_basename.replace('.xlsx','') + output_file_suffix
        num_peptides = int(pool_spot_counts_xlsx_file_basename.split('_')[0].replace('peptides',''))
        if Path(new_elispot_design_xlsx_file).is_file():
            data['file'].append(pool_spot_counts_xlsx_file_basename)
            data['successful'].append(True)
            logger.info('%s was already downloaded' % pool_spot_counts_xlsx_file)
            continue
        if num_peptides == 120:
            design_csv_file = STRANDBERG_TEMPLATES_DIR + '/120peptides_11perpool_3x_strandberg_template_semicolon.csv'
        elif num_peptides == 800:
            design_csv_file = STRANDBERG_TEMPLATES_DIR + '/800peptides_28perpool_3x_strandberg_template_semicolon.csv'
        else:
            raise Exception('Unsupported number of peptides: %i' % num_peptides)
        is_successful = deconvolve_strandberg_heuristic(
            driver=driver,
            pool_spot_counts_xlsx_file=pool_spot_counts_xlsx_file,
            design_csv_file=design_csv_file,
            positivity_criterion=positivity_criterion,
        )
        if is_successful:
            # Rename file
            old_elispot_design_xlsx_file = output_dir + '/elispot_design.xlsx'
            os.rename(old_elispot_design_xlsx_file,  new_elispot_design_xlsx_file)
            logger.info("%s successfully downloaded" % pool_spot_counts_xlsx_file_basename)
        else:
            logger.info("%s failed" % pool_spot_counts_xlsx_file_basename)
        data['file'].append(pool_spot_counts_xlsx_file_basename)
        data['successful'].append(is_successful)
    df_executions = pd.DataFrame(data)
    df_executions.to_csv(output_dir + '/execution_logs.tsv', sep='\t', index=False)


def convert_to_deconvolution_results_helper(output_dir):
    for xlsx_file in glob.glob(output_dir + '/*.xlsx'):
        df_hit_wells = pd.read_excel(xlsx_file, sheet_name='Positive Wells')
        df_positive_antigens = pd.read_excel(xlsx_file, sheet_name='Positive Antigens')
        deconvolution_results_xlsx_file_basename = os.path.basename(xlsx_file)
        deconvolution_result = DeconvolutionResult()
        all_hit_pool_ids = df_hit_wells['Well.nr'].values.tolist()
        for _, row in df_positive_antigens.iterrows():
            peptide_id = 'peptide_%i' % int(row['Antigen'])
            estimated_peptide_spot_count = float(row['Estimate'])
            design_tsv_file = CONFIGURATIONS_DIR + '/' + '_'.join(deconvolution_results_xlsx_file_basename.split('_')[0:7]) + '.tsv'
            peptides_tsv_file = PEPTIDES_DIR + '/' + '_'.join(deconvolution_results_xlsx_file_basename.split('_')[0:3]) + '.tsv'
            df_assignments = pd.read_csv(design_tsv_file, sep='\t')
            df_peptides = pd.read_csv(peptides_tsv_file, sep='\t')
            peptide_sequence = df_peptides.loc[df_peptides['peptide_id'] == peptide_id, 'epitope'].values.tolist()[0]
            pool_ids = df_assignments.loc[df_assignments['peptide_id'] == peptide_id, 'pool_id'].values.tolist()
            hit_pool_ids = []
            for pool_id in pool_ids:
                if pool_id in all_hit_pool_ids:
                    hit_pool_ids.append(pool_id)
            deconvolution_result.add_peptide(
                peptide_id=peptide_id,
                peptide_sequence=peptide_sequence,
                estimated_peptide_spot_count=estimated_peptide_spot_count,
                label=DeconvolutionLabels.CANDIDATE_HIT,
                hit_pool_ids=hit_pool_ids
            )
        output_file = output_dir + '/' + deconvolution_results_xlsx_file_basename.replace('.xlsx', '.tsv')
        deconvolution_result.to_dataframe().to_csv(output_file, sep='\t', index=False)


def deconvolve_pool_spot_counts():
    pool = mp.Pool(processes=NUM_PROCESSES)

    # 300immunogenic_5nonimmunogenic_1dispersion_fnr0.00
    pool.apply_async(
        run_strandberg_heuristic,
        args=[POOL_SPOT_COUNTS_DIR + "/300immunogenic_30nonimmunogenic_1dispersion_fnr0.00/strandberg/upload_ready",
              OUTPUT_DIR + '/300immunogenic_30nonimmunogenic_1dispersion_fnr0.00/strandberg/strandberg_background_subtracted',
              'background_subtracted',
              '_strandberg-background-subtracted.xlsx']
    )
    pool.apply_async(
        run_strandberg_heuristic,
        args=[POOL_SPOT_COUNTS_DIR + "/300immunogenic_30nonimmunogenic_1dispersion_fnr0.00/strandberg/upload_ready",
              OUTPUT_DIR + '/300immunogenic_30nonimmunogenic_1dispersion_fnr0.00/strandberg/strandberg_empirical',
              'empirical',
              '_strandberg-empirical.xlsx']
    )
    pool.apply_async(
        run_strandberg_heuristic,
        args=[POOL_SPOT_COUNTS_DIR + "/300immunogenic_30nonimmunogenic_1dispersion_fnr0.00/strandberg/upload_ready",
              OUTPUT_DIR + '/300immunogenic_30nonimmunogenic_1dispersion_fnr0.00/strandberg/strandberg_bayesian',
              'bayesian',
              '_strandberg-bayesian.xlsx']
    )
    pool.close()
    pool.join()


def convert_to_deconvolution_results():
    pool = mp.Pool(processes=NUM_PROCESSES)

    # 300immunogenic_5nonimmunogenic_30dispersion_fnr0.00
    pool.apply_async(
        convert_to_deconvolution_results_helper,
        args=[OUTPUT_DIR + '/300immunogenic_30nonimmunogenic_1dispersion_fnr0.00/strandberg/strandberg_background_subtracted']
    )
    pool.apply_async(
        convert_to_deconvolution_results_helper,
        args=[OUTPUT_DIR + '/300immunogenic_30nonimmunogenic_1dispersion_fnr0.00/strandberg/strandberg_empirical']
    )
    pool.apply_async(
        convert_to_deconvolution_results_helper,
        args=[OUTPUT_DIR + '/300immunogenic_30nonimmunogenic_1dispersion_fnr0.00/strandberg/strandberg_bayesian']
    )

    pool.close()
    pool.join()


if __name__ == "__main__":
    deconvolve_pool_spot_counts()
    convert_to_deconvolution_results()
