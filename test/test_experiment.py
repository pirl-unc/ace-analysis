import numpy as np
import pandas as pd
from .data import get_data_path
from acesim.solver_ace import AceSolver
from acesim.solver_randomized_design import RandomizedDesignSolver
from acesim.solver_repeated_design import RepeatedDesignSolver
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
    randomized_solver = RandomizedDesignSolver(name='random', random_seed=1)
    repeated_solver = RepeatedDesignSolver(name='repeated')
    
    # Instantiate experiment
    experiment = Experiment(
        experiment_id=1,
        num_peptides=25,
        num_positives=5,
        num_peptides_per_pool=5,
        coverage=3,
        solvers=[ace_solver, randomized_solver, repeated_solver],
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


def test_experiment_2():
    # Instantiate solvers
    ace_solver = AceSolver(
        name='ace',
        cluster_peptides=True,
        random_seed=42,
        mode='golfy',
        trained_model_file=get_data_path(name='seq_sim_trained_model.pt')
    )
    randomized_solver = RandomizedDesignSolver(name='random', random_seed=1)
    repeated_solver = RepeatedDesignSolver(name='repeated')

    # Instantiate experiment
    experiment = Experiment(
        experiment_id=1,
        num_peptides=90,
        num_positives=9,
        num_peptides_per_pool=9,
        coverage=3,
        solvers=[ace_solver, randomized_solver, repeated_solver],
        random_effects=False,
        df_ref_peptides=pd.read_csv(get_data_path(name='iedb_mmer_all.csv')),
        mu_immunogenic=100.0,
        mu_nonimmunogenic=10.0,
        dispersion_factor=2.0,
        method='threshold',
        peptide_sampling_method='alanine_scanning',
        alpha=0.05,
        num_processes=1
    )

    # Run experiment
    experiment.run(10)
