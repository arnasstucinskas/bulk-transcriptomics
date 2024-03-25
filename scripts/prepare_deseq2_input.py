import pandas as pd
import os

# Extract sample names from the input paths
samples = [os.path.basename(path).split('_best_strand.txt')[0] for path in snakemake.input.counts]
conditions_map = {
    "HBR": "HBR",
    "UHRR": "UHRR"
}
conditions = [conditions_map[sample.split('-')[1]] for sample in samples]
count_matrix = pd.DataFrame()

for path, sample in zip(snakemake.input.counts, samples):
    if os.path.exists(path):
        # Read the header to get the column name for counts
        header = pd.read_csv(path, sep='\t', header=0, nrows=0, skiprows=1)
        last_column_name = header.columns[-1]

        # Load counts data
        df = pd.read_csv(path, sep='\t', header=0, usecols=['Geneid', last_column_name], skiprows=1)
        df = df.rename(columns={'Geneid': 'gene_id', last_column_name: sample})
        df = df.set_index('gene_id')
        
        # Concatenate this sample's counts into the count matrix
        count_matrix = pd.concat([count_matrix, df], axis=1)

# Output the count matrix to CSV
count_matrix.to_csv(snakemake.output.count_matrix)

# Prepare and output the column data
col_data = pd.DataFrame({
    'sample': samples,
    'condition': conditions
})

col_data.to_csv(snakemake.output.col_data, index=False)
