"""
The purpose of this python3 script is to implement the PrecomputedConfigurator dataclass
"""


import random
import pandas as pd
from collections import defaultdict
from dataclasses import dataclass, field
from acelib.block_assignment import BlockAssignment
from acelib.types import *
from .configurator import Configurator


@dataclass(frozen=True)
class StrandbergConfigurator(Configurator):

    def generate_assignment(
            self,
            strandberg_template_csv_file: str,
            peptide_tsv_file: str
    ) -> Tuple[BlockAssignment, PeptidePairs]:
        df_template = pd.read_csv(strandberg_template_csv_file)
        df_peptides = pd.read_csv(peptide_tsv_file, sep='\t')
        assignments = {
            1: defaultdict(list),
            2: defaultdict(list),
            3: defaultdict(list)
        }
        plate_ids = {}
        peptides_count = {}
        for _, row in df_template.iterrows():
            if row['Pool'] == 'Ctrl':
                continue
            pool_id = row['Pool']
            well_id = row['Position']
            plate_ids[pool_id] = (1, well_id)
            num_peptides_per_pool = len(df_template.columns.values.tolist()) - 2
            for i in range(0, num_peptides_per_pool):
                if i == 0:
                    col_header = 'Antigen'
                else:
                    col_header = 'Antigen.%i' % i
                if not pd.isna(row[col_header]):
                    peptide_id = 'peptide_%i' % int(row[col_header])
                    peptide_sequence = df_peptides.loc[df_peptides['peptide_id'] == peptide_id, 'epitope'].values.tolist()[0]
                    if peptide_id not in peptides_count.keys():
                        coverage = 1
                        peptides_count[peptide_id] = 1
                    else:
                        coverage = peptides_count[peptide_id] + 1
                        peptides_count[peptide_id] = peptides_count[peptide_id] + 1
                    assignments[coverage][pool_id].append((peptide_id, peptide_sequence))
        return BlockAssignment(assignments=assignments, plate_ids=plate_ids), []
