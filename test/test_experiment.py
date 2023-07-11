import pandas as pd
import pkg_resources
from .data import get_data_path
from acesim.experiment import Experiment
from acesim.solver_ace import AceSolver
from acesim.solver_bogey import BogeySolver


def test_experiment_ace():
    ace_solver = AceSolver(
        name='ace',
        cluster_peptides=True,
        random_seed=1,
        mode='golfy',
        trained_model_file=pkg_resources.resource_filename('acelib', 'resources/models/seq_sim_trained_model.pt')
    )
    bogey_solver = BogeySolver(
        name='bogey'
    )
    experiment = Experiment(
        num_peptides=25,
        num_positives=2,
        num_peptides_per_pool=5,
        coverage=3,
        solvers=[ace_solver, bogey_solver],
        peptide_scan=True,
        df_ref_peptides=pd.read_csv(get_data_path(name='iedb_mmer_all.csv'))
    )
    experiment.run(num_iterations=10)
