"""
The purpose of this python3 script is to aggregate the results TSV files
into one TSV file.
"""


import pandas as pd
import glob


DATA_DIR = '/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/02_benchmark_ace_runtime_memory'
OUTPUT_DIR = '/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/02_benchmark_ace_runtime_memory'


if __name__ == "__main__":
    df_all = pd.DataFrame()
    for csv_file in glob.glob(DATA_DIR + '/*results*.csv'):
        df = pd.read_csv(csv_file)
        df_all = pd.concat([df_all, df])
    df_all.to_csv(OUTPUT_DIR + '/runtime_memory_benchmark_experiment_results_merged.tsv',
                  index=False,
                  sep='\t')
