import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from acesim.solver_ace import AceSolver
#from acesim.solver_bogey import BogeySolver
from acesim.experiment import Experiment



if __name__ == "__main__":
    # Instantiate solvers
    ace_solver = AceSolver(name='ace', cluster_peptides=True, random_seed=42, mode='golfy', trained_model_file='models/ace_model.pt')
    
    # Instantiate experiment
    experiment = Experiment(
    num_peptides=25,
    num_positives=5,
    num_peptides_per_pool=5,
    coverage=3,
    solvers=[],
    random_effects=False,
    df_ref_peptides=pd.read_csv('../refs/iedb_mmer_all.csv'),
    mu_immunogenic=100.0,
    mu_nonimmunogenic=10.0,
    dispersion_factor=2.0,
    method='threshold',
    alpha=0.05,
    num_processes=1
    )

    # Run experiment
    experiment.run(1000)
