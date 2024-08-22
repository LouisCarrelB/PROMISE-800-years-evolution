import pandas as pd

# Define the gene name and file paths
gene_name = "ENSG00000107643"
GENE = "DATA/" + gene_name + "/"

# Load the path_table and s_exon table
path_table = pd.read_csv(GENE + "thoraxe/path_table.csv")
s_exon = pd.read_csv(GENE + "thoraxe/s_exon_table.csv")

# Iterate over each row in the path_table
for index, row in path_table.iterrows():
    gene_id = row['GeneID']
    transcript_id = row['TranscriptIDCluster'].replace('/', '_')  # Replace '/' with '_'
    fasta_filename = f"{gene_id}-{transcript_id}.fasta"
    
    # Create a new FASTA file
    with open(fasta_filename, 'w') as fasta_file:
        # Write the FASTA header
        fasta_file.write(f">{gene_id}-{transcript_id}\n")
        
        # Get the list of exon IDs from the Path column, similar to the build_gene_dict function
        exon_list = row['Path'].split('/')
        exon_list = [exon for exon in exon_list if exon not in ['start', 'stop'] and not exon.startswith('0_')]
        
        # Initialize the sequence to be written to the FASTA file
        final_sequence = ''
        
        # Loop through each exon ID in the list
        for i, exon_id in enumerate(exon_list):
            # Find the corresponding sequence in the s_exon table
            matching_rows = s_exon[s_exon['S_exonID'] == exon_id]
            
            if not matching_rows.empty:
                exon_sequence = matching_rows.iloc[0]['S_exon_Sequence']
                if pd.notna(exon_sequence):
                    final_sequence += exon_sequence
                else:
                    # Print the position and exon ID if sequence is missing
                    print(f"Missing sequence for Exon {exon_id} at position {i+1} in {gene_id}-{transcript_id}")
                    final_sequence += 'NNNN'  # Add placeholder if sequence is missing
            else:
                # Print the position and exon ID if exon ID is not found
                print(f"Exon {exon_id} not found at position {i+1} in {gene_id}-{transcript_id}")
                final_sequence += 'NNNN'  # Add placeholder if exon is not found
        
        # Write the concatenated sequence to the FASTA file
        fasta_file.write(final_sequence + '\n')

# No other prints will be executed