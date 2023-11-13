"""
The purpose of this python3 script is to sample peptides.
"""


import os
import pandas as pd
from acesim.experiment import Experiment


CAMERON_2013_PEPTIDES_CSV_FILE = "/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/raw/sars_cov2_spike_protein_9mers/sars_cov2_spike_protein_9mers_tiling_window.csv"
OUTPUT_DIR = "/datastore/lbcfs/collaborations/pirl/members/jinseok/projects/project_ace/data/processed/04_real_world_datasets/sars_cov2_spike_protein_9mers/01_sampled_peptides"
NUM_REPLICATES = 100


if __name__ == "__main__":
    print("Started")
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    df_ref_peptides = pd.read_csv(CAMERON_2013_PEPTIDES_CSV_FILE)
    for i in range(0, NUM_REPLICATES):
        num_peptides = len(df_ref_peptides)
        num_peptides_immunogenic = len(df_ref_peptides.loc[df_ref_peptides['binding'] == 1,:])
        df_shuffled = df_ref_peptides.sample(frac=1, random_state=Experiment.generate_random_seed()).reset_index(drop=True)
        df_shuffled = df_shuffled.loc[:,['epitope','binding']]
        df_shuffled['peptide_id'] = ['peptide_%i' % (j + 1) for j in range(0, num_peptides)]
        df_shuffled.to_csv('%s/%ipeptides_%iimmunogenic_rep%i.tsv' %
                           (OUTPUT_DIR, num_peptides, num_peptides_immunogenic, i + 1),
                           sep='\t', index=False)
    print("Finished successfully")
