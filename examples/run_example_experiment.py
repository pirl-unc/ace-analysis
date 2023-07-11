import pandas as pd
import multiprocessing as mp
from acesim import Experiment
from dataclasses import dataclass, field


NUM_PROCESSES = 48
NUM_ITERATIONS = 10000
CONFIGURATIONS_TSV_FILE = 'some path'


@dataclass(frozen=True)
class EliSpotConfiguration:
    num_peptides: int
    num_pepties_per_pool: int
    num_coverage: int


    # Define the assay parameters
    num_peptides: int
    num_positives: int
    num_peptides_per_pool: int
    coverage: int
    
    # Define the solver parameters
    solvers: List[Solver]

    # Define the simulation parameters
    random_effects: bool = False
    peptide_scan: bool = False
    mu_immunogenic: float = 100.0
    mu_nonimmunogenic: float = 10.0
    dispersion_factor: float = 1.0

    # Define the hardware parameters
    num_processes: int = 1
    random_state: int = Experiment.generate_random_seed()

    # Define the reference database
    _df_ref_peptides: pd.DataFrame = None


def worker(experiment: Experiment) -> pd.DataFrame:
    """
    Performs one experiment using Experiment class.
    """
    # Step 1. Create an Experiment object.
    experiment = Experiment(
        num_peptides=elispot_configuration.num_peptides,
        num_peptides_per_pool=elispot_configuration.
    )

    # Step 1. Sample peptides
    df_peptides = experiment.sample_peptides()

    # Step 2. Run ACE to generate an ELISpot configuration
    
    
    # Step 3. Run Bogey

    # Step 4. Simulate spot counts for each peptide
    peptide_spot_counts = experiment.simulate_peptide_spot_counts()

    # Step 5. Simulate spot counts for pools
    pool_spot_counts = experiment.simulate_pool_spot_counts()

    # Step 6. Identify hit pool IDs
    hit_pool_ids = blah



if __name__ == "__main__":
    df_experiment_configs = pd.read_csv(CONFIGURATIONS_TSV_FILE, sep='\t')
    for index, row in df_experiment_configs.iterrows():
        num_peptides = int(row['num_petides'])
        num_peptides_per_pool = int(row['num_peptides_per_pool'])
        num_coverage = int(row['num_coverage'])