"""
The purpose of this python3 script is to aggregate the results TSV files
into one TSV file.
"""


import pandas as pd
import glob


DATA_DIR = '/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/01_benchmark_ace/real_world_datasets/tan_et_al_inf_gen_evol_2021'
OUTPUT_DIR = '/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/01_benchmark_ace/real_world_datasets/tan_et_al_inf_gen_evol_2021'


if __name__ == "__main__":
    df_all = pd.DataFrame()
    for tsv_file in glob.glob(DATA_DIR + '/*/*results*.tsv'):
        print(tsv_file)
        df = pd.read_csv(tsv_file, sep='\t')
        df.drop(columns=['preferred_peptide_pairs', 'positive_peptide_sequences'], inplace=True)
        df_all = pd.concat([df_all, df])
    df_all.to_csv(OUTPUT_DIR + '/tan2021_dataset_experiment_results_merged.tsv',
                  index=False,
                  sep='\t')
