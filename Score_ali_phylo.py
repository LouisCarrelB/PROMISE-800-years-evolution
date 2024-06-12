import pandas as pd
import json
import subprocess
import os

def read_msa_fasta(file_path):
    """
    Read a MSA FASTA file and extract the ID, pA, and pB from each sequence.
    Args:
        file_path (str): Path to the FASTA file.
    Returns:
        DataFrame: Contains columns ID, pA, and pB.
    """
    sequences = []
    try:
        with open(file_path, 'r') as file:
            for line in file:
                if line.startswith('>'):
                    header = line.strip()
                    id_part = header.split(' ')[0][1:]
                    try:
                        pA = float(header.split('pA=')[1].split(' ')[0])
                        pB = float(header.split('pB=')[1].split(' ')[0])
                        sequences.append([id_part, pA, pB])
                    except IndexError:
                        print(f"Index error in sequence {id_part}, header: {header}")
    except FileNotFoundError:
        print(f"File not found: {file_path}")
        return pd.DataFrame(columns=['ID', 'pA', 'pB'])
    
    return pd.DataFrame(sequences, columns=['ID', 'pA', 'pB'])

def merge_msas(file1, file2):
    """
    Merge data from two MSA files into a single DataFrame.
    Args:
        file1, file2 (str): Paths to the FASTA files.
    Returns:
        DataFrame: Merged data with columns ID, pA, and pB.
    """
    df1 = read_msa_fasta(file1)
    df2 = read_msa_fasta(file2)
    return pd.concat([df1, df2], ignore_index=True)

def format_species_name(species_name):
    """
    Format the species name to start with a capital letter and replace spaces with underscores.
    The first word starts with an uppercase letter and the second word starts with a lowercase letter.
    Args:
        species_name (str): The original species name.
    Returns:
        str: The formatted species name.
    """
    parts = species_name.split()
    if len(parts) == 2:
        return f"{parts[0].capitalize()}_{parts[1].lower()}"
    return species_name.capitalize()

def add_species_column(df, species_mapping):
    """
    Add a 'Species' column to a DataFrame using a mapping dictionary and format the names.
    Args:
        df (DataFrame): Original DataFrame containing sequence IDs.
        species_mapping (dict): Dictionary mapping IDs to species names.
    Returns:
        DataFrame: Modified DataFrame with 'Species' column added.
    """
    df['Species'] = df['ID'].map(species_mapping).fillna('Unknown')
    df['Species'] = df['Species'].apply(format_species_name)
    return df

def execute_pastml(pastml_command):
    """
    Execute a PasteML command using subprocess.
    Args:
        pastml_command (str): Command for running PasteML.
    """
    try:
        subprocess.call(pastml_command, shell=True)
    except Exception as e:
        print(f"Error executing PasteML command: {e}")

def create_species_mapping(csv_path):
    """
    Create a unique species mapping dictionary from a CSV file.
    Args:
        csv_path (str): Path to the CSV file containing gene and species information.
    Returns:
        dict: Dictionary mapping unique gene identifiers to species names.
    """
    try:
        # Load data from CSV
        data = pd.read_csv(csv_path, usecols=['GeneID', 'Species'])
        
        # Drop duplicates to ensure unique combinations
        unique_data = data.drop_duplicates(subset=['GeneID', 'Species'])
        
        # Create dictionary mapping geneID to Species
        species_mapping = pd.Series(unique_data.Species.values, index=unique_data.GeneID).to_dict()
        
        return species_mapping
    except FileNotFoundError:
        print(f"File not found: {csv_path}")
        return {}
    except Exception as e:
        print(f"Error reading or processing CSV file: {e}")
        return {}

def main():
    base_path = '/Users/louiscarrel/Documents/Alignement_Project/largescale_kinase/DATA/ENSG00000107643'
    msa_files = [
        os.path.join(base_path, 'New_alignement/msa_s_exon_17_0.fasta'),
        os.path.join(base_path, 'inter/Alt_17_0.fasta')
    ]
    result_df = merge_msas(*msa_files)

    species_csv_path = os.path.join(base_path, 'thoraxe/s_exon_table_a3m.csv')
    species_dict = create_species_mapping(species_csv_path)
    result_df = add_species_column(result_df, species_dict)
    
    result_df = result_df[result_df['Species'] != 'Unknown']
    result_df.reset_index(drop=True, inplace=True)
    result_df = result_df[['pA', 'Species']].rename(columns={'pA': 'Index'})

    output_tsv_path = os.path.join(base_path, 'result_with_species_with_proba.tsv')
    result_df.to_csv(output_tsv_path, index=False, sep='\t')

    # pastml_command = f"pastml --tree {os.path.join(base_path, 'Species_list_time_tree.nwk')} --data {output_tsv_path} --html_compressed {os.path.join(base_path, 'map.html')}"
    # execute_pastml(pastml_command)

if __name__ == "__main__":
    main()
