"""
The purpose of this python3 script is to implement utility functions
"""


from golfy.design import Design
from golfy.types import SpotCounts
from acelib.block_assignment import BlockAssignment
from acelib.block_design import BlockDesign
from typing import Dict, Iterable, Mapping, Tuple


PeptideIndicesToPeptideIds = Mapping[int, str]


def convert_block_assignment_to_golfy_design(
        block_assignment: BlockAssignment
) -> Tuple[Design, PeptideIndicesToPeptideIds]:
    """
    Converts a BlockAssignment object to golfy Design object.

    Parameters
    ----------
    block_assignment        :   BlockAssignment object.

    Returns
    -------
    design                  :   Golfy Design object.
    peptide_idx_to_id_dict  :   Dictionary where
                                key = peptide index
                                value = peptide ID
    """
    # Step 1. Create dictionaries of peptide indices to IDs (and vice and versa)
    peptide_idx_to_id_dict = {}
    peptide_id_to_idx_dict = {}
    idx = 0
    for peptide_id in block_assignment.peptide_ids:
        peptide_idx_to_id_dict[idx] = peptide_id
        peptide_id_to_idx_dict[peptide_id] = idx
        idx += 1

    # Step 2. Create assignments for golfy Design class
    assignments = {}
    num_peptides_per_pool = -1
    for coverage in block_assignment.assignments.keys():
        assignments[coverage-1] = {}
        for pool in block_assignment.assignments[coverage].keys():
            assignments[coverage-1][pool-1] = []
            for peptide_id, peptide_sequence in block_assignment.assignments[coverage][pool]:
                assignments[coverage-1][pool-1].append(peptide_id_to_idx_dict[peptide_id])
            if num_peptides_per_pool < len(block_assignment.assignments[coverage][pool]):
                num_peptides_per_pool = len(block_assignment.assignments[coverage][pool])

    # Step 3. Create a golfy Design object
    design = Design(
        num_peptides=len(block_assignment.peptide_ids),
        max_peptides_per_pool=num_peptides_per_pool,
        num_replicates=len(block_assignment.assignments.keys()),
        allow_extra_pools=False,
        invalid_neighbors=[],
        preferred_neighbors=[],
        assignments=assignments
    )
    return design, peptide_idx_to_id_dict


def convert_spot_counts_to_golfy_spotcounts(
        spot_counts: Dict[int, int], 
        block_assignment: BlockAssignment
) -> SpotCounts:
    """
    Converts spot counts to golfy SpotCounts.

    Parameters
    ----------
    spot_counts             :   Dictionary where
                                key = pool ID
                                value = spot count
    block_assignment        :   BlockAssignment object.

    Returns
    -------
    counts                  :   golfy SpotCounts.
    """


    counts = {}
    for coverage in block_assignment.assignments.keys():
        counts[coverage-1] = {}
        for pool in block_assignment.assignments[coverage].keys():
            counts[coverage-1][pool-1] = spot_counts[pool]
    return counts
