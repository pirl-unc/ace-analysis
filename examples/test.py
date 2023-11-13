from acelib.main import run_ace_deconvolve
from acelib.block_assignment import BlockAssignment
import pandas as pd

POOL_SPOT_COUNTS_FILE = "/Users/leework/Documents/Research/projects/project_ace/data/processed/test/90peptides_1immunogenic_rep23_9perpool_3x_ace-s_design_98625556_pool_spot_counts.tsv"
CONFIG_FILE = "/Users/leework/Documents/Research/projects/project_ace/data/processed/test/90peptides_1immunogenic_rep23_9perpool_3x_ace-s_design_98625556.tsv"

d = run_ace_deconvolve(
    df_readout=pd.read_csv(POOL_SPOT_COUNTS_FILE, sep='\t'),
    block_assignment=BlockAssignment.load_from_dataframe(df_assignments=pd.read_csv(CONFIG_FILE, sep='\t')),
    method='cem',
    min_coverage=3,
    min_pool_spot_count=100
)
df = d.to_dataframe()
print(df.loc[df['deconvolution_result'] != 'not_a_hit',:])
