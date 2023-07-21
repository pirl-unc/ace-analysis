import numpy as np
import pandas as pd
from .data import get_data_path
from acesim.solver_ace import AceSolver
from acesim.experiment import Experiment


def test_experiment_1():
    # Instantiate solvers
    ace_solver = AceSolver(
        name='ace',
        cluster_peptides=True,
        random_seed=42,
        mode='golfy',
        trained_model_file=get_data_path(name='seq_sim_trained_model.pt')
    )
    
    # Instantiate experiment
    experiment = Experiment(
        experiment_id=1,
        num_peptides=25,
        num_positives=5,
        num_peptides_per_pool=5,
        coverage=3,
        solvers=[ace_solver],
        random_effects=False,
        df_ref_peptides=pd.read_csv(get_data_path(name='iedb_mmer_all.csv')),
        mu_immunogenic=100.0,
        mu_nonimmunogenic=10.0,
        dispersion_factor=2.0,
        method='threshold',
        alpha=0.05,
        num_processes=1
    )

    # Run experiment
    experiment.run(10)
