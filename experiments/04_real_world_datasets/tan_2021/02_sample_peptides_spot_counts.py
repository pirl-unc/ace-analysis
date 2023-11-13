"""
The purpose of this python3 script is to sample spot counts for peptides.
"""


import glob
import os
import pandas as pd
from acesim.experiment import Experiment


SAMPLED_PEPTIDES_DIR = "/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/04_real_world_datasets/tan_2021/01_sampled_peptides"
OUTPUT_DIR = "/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/04_real_world_datasets/tan_2021/02_sampled_peptides_spot_counts"
MU_IMMUNOGENIC = 300.0
MU_NONIMMUNOGENIC = 30.0
DISPERSION_FACTOR = 1.0


if __name__ == "__main__":
    print("Started")
    output_dir_ = OUTPUT_DIR + '/300immunogenic_30nonimmunogenic_1dispersion'
    if not os.path.exists(output_dir_):
        os.makedirs(output_dir_)
    for tsv_file in glob.glob(SAMPLED_PEPTIDES_DIR + '/*.tsv'):
        tsv_file_basename = os.path.basename(tsv_file)
        df_peptides = pd.read_csv(tsv_file, sep='\t')
        df_peptides_spot_counts = Experiment.simulate_peptide_spot_counts(
            df_peptides=df_peptides,
            num_coverage=3,
            random_effects=True,
            mu_immunogenic=MU_IMMUNOGENIC,
            mu_nonimmunogenic=MU_NONIMMUNOGENIC,
            dispersion_factor=DISPERSION_FACTOR
        )
        output_file = output_dir_ + '/' + tsv_file_basename.replace('.tsv', '') + '_3x_spot_counts.tsv'
        df_peptides_spot_counts.to_csv(output_file, sep='\t', index=False)
    print("Finished successfully")
