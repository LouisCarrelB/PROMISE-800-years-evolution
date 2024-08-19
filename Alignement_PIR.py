'''
This script takes as input the name of a gene and compares thorax transcripts with protein sequences in the a3m. The a3m must be in the folder named after the gene, which contains the "thorax" folder, 
otherwise, this needs to be modified in the script. The output is a CSV file named a3m_to_PIR, which contains the protein sequences from the a3m and the edit distances with the closest transcripts.
We also take two other arguments: the second is the number of sequences examined in the a3m (put 1 if you want only the query sequence to be analyzed), and the third is the maximum edit distance for the matches to be included in the output.
Author: Carrel-Billiard Louis
Date: 03/07/2024
'''

import editdistance
import re
import glob
from Bio import SeqIO, AlignIO
import requests
import csv
import pandas as pd
import sys

def read_pir_file(file_path):
    '''
    Reads a PIR file and returns the lines.

    Parameters:
    file_path (str): Path to the PIR file.

    Returns:
    list: List of lines in the file.
    '''
    with open(file_path, 'r') as file:
        lines = file.readlines()
    return lines

def find_sequence_in_pir(sequence, pir_lines, max_distance):
    '''
    Finds a sequence in a PIR file with a maximum edit distance.

    Parameters:
    sequence (str): The target sequence.
    pir_lines (list): List of lines from the PIR file.
    max_distance (int): Maximum edit distance to consider a matching sequence.

    Returns:
    tuple: Gene identifiers, transcript identifiers, and found sequences.
    '''
    first_identifier = []
    second_identifier = []
    found_sequences = []

    for line in pir_lines:
        if line.startswith('>'):
            current_sequence = line.strip()[1:]
        else:
            distance = editdistance.eval(sequence, line.strip())
            if distance <= max_distance:
                found_sequences.append((current_sequence, distance))
                # Extract the part between the first and next space
                match = re.search(r';([^ ]+)\s([^ ]+)', current_sequence)
                if match:
                    first_identifier.append(match.group(1))
                    second_identifier.append(match.group(2))
    return first_identifier, second_identifier, found_sequences

def extract_uni_id(sequence_id):
    '''
    Extracts the UniProt identifier from the sequence identifier.

    Parameters:
    sequence_id (str): Sequence identifier.

    Returns:
    str: UniProt identifier.
    '''
    pattern = re.compile(r"_(.+)")
    match = pattern.search(sequence_id)
    if match:
        return match.group(1)
    else:
        return sequence_id

def get_organism_name(uniprot_id):
    '''
    Gets the organism name from the UniProt identifier.

    Parameters:
    uniprot_id (str): UniProt identifier.

    Returns:
    str: Scientific name of the organism.
    '''
    url = f"https://www.ebi.ac.uk/proteins/api/proteins/{uniprot_id}"
    response = requests.get(url)

    if response.status_code == 200:
        protein_data = response.json()
        organism_names = protein_data.get("organism", {}).get("names", [])

        for name_entry in organism_names:
            if name_entry.get("type", "") == "scientific":
                return name_entry.get("value", "")

    return None

def create_gene_species_dict(file_path):
    '''
    Creates a dictionary associating GeneIDs with species from a CSV file.

    Parameters:
    file_path (str): Path to the CSV file.

    Returns:
    dict: Dictionary linking GeneIDs to species.
    '''
    s_exons_df = pd.read_csv(file_path, delimiter=',')
    unique_gene_ids = s_exons_df['GeneID'].astype(str).unique()
    gene_species_dict = {}
    for gene_id in unique_gene_ids:
        gene_species_dict[gene_id] = list(s_exons_df[s_exons_df['GeneID'] == gene_id]['Species'].unique())
    return gene_species_dict

if __name__ == "__main__":
    '''
    Main entry point of the script. Reads command line arguments,
    searches for sequences in a PIR file, and generates a CSV file with the results.
    '''
    csv_data = []
    if len(sys.argv) < 2:
        print("Usage: python Alignement_PIR.py <gene_name> <nbr_of_sequence_to_look_in_a3m> <max_edit_distance>")
        sys.exit(1)
        
    gene_name = sys.argv[1]
    NUM_a3m = int(sys.argv[2])
    max_distance = int(sys.argv[3])



    ############################### Can be modify following the config of you folder 
    GENE = "DATA/" + gene_name + "/"
    ###############################



    pir_file_path = GENE + 'thoraxe/phylosofs/transcripts.pir'  
    s_exons = GENE + "thoraxe/s_exon_table.csv"

    #here the a3m must be in your "GENE" folder
    for a3m_file in glob.glob(GENE + "*.a3m"):
        sequences_a3m = list(SeqIO.parse(a3m_file, "fasta"))[:NUM_a3m]
        for sequence_a3m in sequences_a3m:
            pattern = re.compile(r"_(.+)")

            target_sequence = sequence_a3m.seq
            sequence_id = sequence_a3m.id

            print(extract_uni_id(sequence_id))
            organism_name = get_organism_name(extract_uni_id(sequence_id))
            pir_lines = read_pir_file(pir_file_path)
            results_ID_g, results_ID_t, results_tot = find_sequence_in_pir(target_sequence, pir_lines, max_distance)
            print("For the sequence:", sequence_id)
            print("Organism:", organism_name)
            if results_ID_g:
                print("Found Sequences:")
                for sequence, distance in results_tot:
                    print(f"{sequence} (Edit distance: {distance})")
                print("Transcript IDs found:")
                for identifier in results_ID_t:
                    print(identifier)
                print("Gene IDs found:")
                for identifier in results_ID_g:
                    print(identifier)
                for i, (sequence, distance) in enumerate(results_tot):
                    transcript_id = results_ID_t[i] if i < len(results_ID_t) else None
                    gene_id = results_ID_g[i] if i < len(results_ID_g) else None

                    csv_data.append({
                        "Sequence ID": sequence_id,
                        "Organism_a3m": organism_name,
                        "Found Sequences": sequence,
                        "Distance": distance,
                        "Transcript IDs": transcript_id,
                        "Gene IDs": gene_id,
                    })
            else:
                print("No corresponding identifiers found.") 

    dic_ensembl = create_gene_species_dict(s_exons)
    csv_data_df = pd.DataFrame(csv_data)
    csv_data_df['organism_pir'] = csv_data_df['Gene IDs'].map(dic_ensembl)
    print(csv_data_df)
    csv_data_df.to_csv(GENE + 'a3m_to_PIR.csv', index=False)
