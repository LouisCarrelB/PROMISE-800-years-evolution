"""
Script : Fasta_for_a3M.py
Description :
    Ce script génère des fichiers FASTA à partir des données de transcrits pour un gène donné.
    Il est utilisé avant l'exécution de ColabFold pour préparer les alignements multiples au
    format A3M. Le script lit deux fichiers CSV : 'path_table.csv' et 's_exon_table.csv',
    qui contiennent des informations sur les chemins d'exons et les séquences d'exons
    respectivement. Il génère ensuite des fichiers FASTA, où chaque fichier correspond à un
    transcrit et contient la séquence concaténée des exons.

    Fonctionnement :
    - Pour chaque gène et transcrit, un dossier est créé dans le répertoire 'results'.
    - Le script itère sur les exons d'un transcrit, concatène leurs séquences, et écrit le
      tout dans un fichier FASTA.
    - En cas de séquence manquante ou d'exon introuvable, le script ajoute une séquence
      placeholder ('NNNN') et signale le problème.

Utilisation :
    python Fasta_for_a3m.py <gene_name>

    Où <gene_name> correspond au nom du gène pour lequel les FASTA doivent être générés.

Exemple :
    python3 Fasta_for_a3m.py ENSG00000107643
"""



import os
import sys
import pandas as pd

def generate_fasta_for_transcripts(gene_name):
    # Define the file paths
    GENE = f"DATA/{gene_name}/"

    # Load the path_table and s_exon table
    path_table = pd.read_csv(GENE + "thoraxe/path_table.csv")
    s_exon = pd.read_csv(GENE + "thoraxe/s_exon_table.csv")

    # Iterate over each row in the path_table
    for index, row in path_table.iterrows():
        gene_id = row['GeneID']
        transcript_id = row['TranscriptIDCluster'].replace('/', '_')  # Replace '/' with '_'
        
        # Create directories for the gene and transcript
        gene_dir = os.path.join("results", gene_id)
        transcript_dir = os.path.join(gene_dir, transcript_id)
        os.makedirs(transcript_dir, exist_ok=True)
        
        fasta_filename = os.path.join(transcript_dir, 'transcript.fasta')
        
        # Create a new FASTA file
        with open(fasta_filename, 'w') as fasta_file:
            # Write the FASTA header
            fasta_file.write(f">{gene_id}-{transcript_id}\n")
            
            # Get the list of exon IDs from the Path column
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

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 Fasta_for_a3m.py <gene_name>")
        sys.exit(1)
    
    gene_name = sys.argv[1]
    generate_fasta_for_transcripts(gene_name)