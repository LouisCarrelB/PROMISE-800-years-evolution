import pandas as pd
import json
import subprocess
import os

def read_msa_fasta(file_path):
    sequences = []
    try:
        with open(file_path, 'r') as file:
            for line in file:
                if line.startswith('>'):
                    header = line.strip()
                    id_part = header.split(' ')[0][1:]
                    try:
                        # Modification ici pour adapter à votre format spécifique
                        pA = float(header.split('pA:')[1].split(',')[0].strip())
                        pB = float(header.split('pB:')[1].split(',')[0].strip())
                        sequences.append([id_part, pA, pB])
                    except IndexError:
                        print(f"Index error in sequence {id_part}, header: {header}")
    except FileNotFoundError:
        print(f"File not Foun:id_part, pA, pB: {file_path}")
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
    




def save_species_list(df, output_path):
    """
    Save a list of all unique species in the DataFrame to a TSV file.
    Args:
        df (DataFrame): DataFrame containing a 'Species' column.
        output_path (str): Path to save the TSV file.
    """
    unique_species = df['Species'].unique()
    species_df = pd.DataFrame(unique_species, columns=['Species'])
    species_df.to_csv(output_path, index=False, sep='\t')
    print(f"Species list saved to {output_path}")


def classify_species_origin(df, file1, file2, output_path):
    """
    Classify species as originating from 'Alt', 'Can', or 'Both' based on sequence IDs from two FASTA files.
    Args:
        df (DataFrame): DataFrame containing columns 'ID' and 'Species'.
        file1, file2 (str): Paths to the FASTA files used to classify origin.
        output_path (str): Path to save the TSV file.
    """
    # Read sequences from each file and extract IDs
    ids1 = set(read_msa_fasta(file1)['ID'])
    ids2 = set(read_msa_fasta(file2)['ID'])
    
    # Check if species in df are in both sets
    species_in_ids1 = set(df[df['ID'].isin(ids1)]['Species'])
    species_in_ids2 = set(df[df['ID'].isin(ids2)]['Species'])
    
    # Intersection gives us species present in both
    species_in_both = species_in_ids1.intersection(species_in_ids2)

    # Function to determine origin
    def determine_origin(species):
        if species in species_in_both:
            return 'Both'
        elif species in species_in_ids1:
            return 'Alt'
        elif species in species_in_ids2:
            return 'Can'
        return 'Unknown'

    # Map species to their origin
    df['Index'] = df['Species'].apply(determine_origin)
    
    # Drop duplicates to keep only unique combinations
    output_df = df[['Species', 'Index']].drop_duplicates().reset_index(drop=True)
    output_df.to_csv(output_path, index=False, sep='\t')
    print(f"File with species origins saved to {output_path}")






def main():
    base_path = '/Users/louiscarrel/Documents/Alignement_Project/largescale_kinase/DATA/ENSG00000107643/Analyze_logo/aligned_sequences/'
    msa_files = [
        os.path.join(base_path, 'can_ali.fasta'),
        os.path.join(base_path, 'alt_ali.fasta')
    ]
    result_df = merge_msas(*msa_files)

    species_csv_path = os.path.join('/Users/louiscarrel/Documents/Alignement_Project/largescale_kinase/DATA/ENSG00000107643/thoraxe/s_exon_table_a3m.csv')
    species_dict = create_species_mapping(species_csv_path)
    result_df = add_species_column(result_df, species_dict)
    
    result_df = result_df[result_df['Species'] != 'Unknown']
    result_df.reset_index(drop=True, inplace=True)

    # Changement ici: garder 'ID' pour l'utilisation dans classify_species_origin
    result_df = result_df[['ID', 'pA', 'Species']]
    result_df.rename(columns={'pA': 'Index'}, inplace=True)

    output_tsv_path = os.path.join(base_path, 'result_with_species_with_proba.tsv')
    result_df.to_csv(output_tsv_path, index=False, sep='\t')

    # Save all species to TSV
    species_list_output_path = os.path.join(base_path, 'all_species_list.tsv')
    save_species_list(result_df, species_list_output_path)

    # Save species origin classification
    species_origin_output_path = os.path.join(base_path, 'species_origin_classification.tsv')
    classify_species_origin(result_df, msa_files[0], msa_files[1], species_origin_output_path)

    # Optionally run PasteML
    # pastml_command = f"pastml --tree {os.path.join(base_path, 'Species_list_time_tree.nwk')} --data {output_tsv_path} --html_compressed {os.path.join(base_path, 'map.html')}"
    # execute_pastml(pastml_command)

if __name__ == "__main__":
    main()
