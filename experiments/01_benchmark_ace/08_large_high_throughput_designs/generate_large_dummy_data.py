import random
import pandas as pd


OUTPUT_DIR = "/Users/leework/Documents/Research/projects/project_ace/data/raw"
AMINO_ACIDS = ['A','R','N','D','C','E','Q','G','H','I',
               'L','K','M','F','P','S','T','W','Y','V']


if __name__ == "__main__":
    data = {
        'Epitope': [],
        'Allele': [],
        'Binding': []
    }

    # Step 1. Append positive sequences
    for _ in range(0, 4000):
        epitope = ''.join(random.choices(AMINO_ACIDS, k=9))
        data['Epitope'].append(epitope)
        data['Allele'].append('HLA-A*02:01')
        data['Binding'].append(1)

    # Step 2. Append negative sequences
    for _ in range(0, 36000):
        epitope = ''.join(random.choices(AMINO_ACIDS, k=9))
        data['Epitope'].append(epitope)
        data['Allele'].append('HLA-A*02:01')
        data['Binding'].append(0)

    df = pd.DataFrame(data)
    df.to_csv(OUTPUT_DIR + '/large_dummy_sequence_dataset.csv', index=False)
