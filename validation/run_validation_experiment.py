"""
The purpose of this python3 script is to 
run an in silico validation experiment for ACE.
"""


import argparse
import random
from baseline_methods import *


def parse_args():
    """
    Initializes the input argument parser.
    """
    arg_parser = argparse.ArgumentParser(
        description="Run an in silico validation experiment for ACE."
    )
    arg_parser.add_argument(
        '--num_peptides',
        dest='num_peptides',
        type=int,
        required=True,
        help="Number of total peptides."
    )
    arg_parser.add_argument(
        '--num_true_positive_peptides',
        dest='num_true_positive_peptides',
        type=int,
        required=True,
        help="Number of true positive peptides."
    )
    arg_parser.add_argument(
        '--num_peptides_per_pool',
        dest='num_peptides_per_pool',
        type=int,
        required=True,
        help="Number of peptides per pool."
    )
    arg_parser.add_argument(
        '--num_coverage',
        dest='num_coverage',
        type=int,
        required=True,
        help="Number of coverage (i.e. repeats per peptide)."
    )
    arg_parser.add_argument(
        '--false_positive_rate',
        dest='false_positive_rate',
        type=float,
        default=0.05,
        required=True,
        help="False positive rate."
    )
    args = arg_parser.parse_args()
    return args
    


"""
Run ACE Using Arguments Passed
"""

def run_ace()
    pass


def identify_hit_peptide_ids(df_config: pd.DataFrame, 
                             hit_pool_ids: dict) -> List[str]:
    """
    Identifies hit peptide IDs.

    Parameters
    ----------
    df_config       :   DataFrame with the following columns:
                        'peptide_id'
                        'pool_id'
                        'coverage_id'
    hit_pool_ids    :   Dictionary where the key is 'pool_id'
                        and the value is 'is_hit'

    Returns
    -------
    hit_peptide_ids :   Hit peptide IDs. 
    """
    
    

if __name__ == "__main__":
    # Step 1. Parse args
    args = parse_args()

    # Step 2. Generate peptide IDs
    peptide_ids = ['peptide_' + str(i) for i in range(1, args.num_peptides + 1)]

    # Step 3. Randomly sample true positive (TP) hit peptide IDs
    peptide_ids_tp = random.sample(peptide_ids, k=args.num_true_positive_peptides)

    # Step 4. Generate ACE configuration
    df_config_ace = run_ace() 

    # Step 5. Generate benchmark configurations
    # Random assignment
    df_config_random_1 = generate_random_elispot_configuration(
        peptide_ids=peptide_ids,
        peptides_per_pool=args.num_peptides_per_pool,
        num_coverage=args.num_coverage
    )

    # Step 6. Perform experiment
    # Generate True Positive Data and Apply Noise
    alpha = 0.05 # False Positive Rate
   

    def generate_empirical_results(df_config, peptide_ids_tp, p=alpha):
        """
        This function needs to take in a given config and true
        peptide labels and output the empirical pool hit observations

        Parameters
        ----------
        df_config       :   DataFrame with the following columns:
                            'peptide_id'
                            'pool_id'
                            'coverage_id'
        peptide_ids_tp  :   List of true positive hit peptide IDs.
        p               :   Probability of success from 1 Bernoulli trial.
        
        Returns
        -------
        Dictionary where the key is 'pool_id' and the value is 'is_hit'
        """
        pools = df_config['pool_id'] #get all the pool IDs
        n_pools = len(pools)
        
        hits_dict = dict(zip(pools, np.zeros(n_pools))) #instantiate the hit pools as zero

        #Get the actual hit pools based on ground truths
        true_pool_hits = [pools[idx] for idx, element in enumerate(df_config['peptide_id']) if element in peptide_ids_tp]
        
        for pool_hit in true_pool_hits:
            hits_dict[pool_hit] = 1
        
        perturb = np.random.binomial(n_pools, p)
        
        for i in range(n_pools):
            if perturb[i] == 1:
                hits_dict[pools[i]] : int(not pools[i])
            else:
                continue
        
        return hits_dict
                
