import pandas as pd
from .data import get_data_path
from acelib.block_assignment import BlockAssignment
from acesim.solver_ace import AceSolver
from acesim.solver_precomputed_design import PrecomputedSolver
from acesim.solver_randomized_design import RandomizedDesignSolver
from acesim.solver_repeated_design import RepeatedDesignSolver
from acesim.experiment import Experiment


def test_experiment_1():
    # Step 1. Instantiate solvers
    ace_solver = AceSolver(
        name='ace',
        cluster_peptides=True,
        mode='golfy',
        trained_model_file=get_data_path(name='trained_model_w_negatives.pt')
    )
    randomized_solver = RandomizedDesignSolver(name='random')
    repeated_solver = RepeatedDesignSolver(name='repeated')
    
    # Step 2. Instantiate experiment
    experiment = Experiment(
        experiment_id=1,
        num_peptides=25,
        num_positives=5,
        num_peptides_per_pool=5,
        coverage=3,
        solvers=[ace_solver, randomized_solver, repeated_solver],
        random_effects=True,
        df_ref_peptides=pd.read_csv(get_data_path(name='held_out_data_w_negatives.csv')),
        mu_immunogenic=100.0,
        mu_nonimmunogenic=10.0,
        dispersion_factor=2.0,
        method='threshold',
        alpha=0.05,
        num_processes=1
    )

    # Step 3. Run experiment
    experiment.run(1)


def test_experiment_2():
    # Step 1. Instantiate solvers
    ace_solver = AceSolver(
        name='ace',
        cluster_peptides=True,
        mode='golfy',
        trained_model_file=get_data_path(name='trained_model_w_negatives.pt')
    )
    randomized_solver = RandomizedDesignSolver(name='random')
    repeated_solver = RepeatedDesignSolver(name='repeated')

    # Step 2. Instantiate experiment
    experiment = Experiment(
        experiment_id=1,
        num_peptides=90,
        num_positives=9,
        num_peptides_per_pool=9,
        coverage=3,
        solvers=[ace_solver, randomized_solver, repeated_solver],
        random_effects=True,
        df_ref_peptides=pd.read_csv(get_data_path(name='held_out_data_w_negatives.csv')),
        mu_immunogenic=100.0,
        mu_nonimmunogenic=10.0,
        dispersion_factor=2.0,
        method='threshold',
        peptide_sampling_method='alanine_scanning',
        alpha=0.05,
        num_processes=1
    )

    # Step 3. Run experiment
    experiment.run(1)


def test_experiment_3():
    # Step 1. Instantiate solvers
    block_assignment = BlockAssignment.read_excel_file(
        excel_file=get_data_path('100peptides_10perpool_3x.xlsx'),
        sheet_name='block_assignment'
    )
    precomputed_solver = PrecomputedSolver(
        name='precomputed',
        block_assignment=block_assignment
    )

    # Step 2. Instantiate experiment
    experiment = Experiment(
        experiment_id=1,
        num_peptides=100,
        num_positives=10,
        num_peptides_per_pool=10,
        coverage=3,
        solvers=[precomputed_solver],
        random_effects=True,
        df_ref_peptides=pd.read_csv(get_data_path(name='held_out_data_w_negatives.csv')),
        mu_immunogenic=100.0,
        mu_nonimmunogenic=10.0,
        dispersion_factor=2.0,
        method='threshold',
        peptide_sampling_method='',
        alpha=0.05,
        num_processes=1
    )

    # Step 3. Run experiment
    experiment.run(1)
