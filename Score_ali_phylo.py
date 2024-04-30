'''
Permet de construire un arbre (avec time tree) et de le mapper avec pasteML en fonction des pA et pB sur les exons alt et can 
'''
import pandas as pd 
import json

def read_msa_fasta(file_path):
    """
    Cette fonction lit un fichier MSA FASTA et extrait les ID, pA, et pB de chaque séquence.

    :param file_path: Chemin vers le fichier FASTA
    :return: DataFrame avec les colonnes ID, pA, et pB
    """
    sequences = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('>'):  # Les lignes d'en-tête commencent par '>'
                header = line.strip()
                # Extraction de l'ID (tout avant le premier espace)
                id_part = header.split(' ')[0][1:]  # Supprime '>'
                
                # Vérifier et extraire pA et pB s'ils existent
                try:
                    pA_part = header.split('pA=')[1].split(' ')[0]
                    pB_part = header.split('pB=')[1].split(' ')[0]
                    pA = float(pA_part)
                    pB = float(pB_part)
                    sequences.append([id_part, pA, pB])
                except IndexError:
                    print(f"Erreur d'index pour la séquence {id_part}, en-tête : {header}")
                    continue  # Continue avec la prochaine ligne si l'en-tête est mal formé

    return pd.DataFrame(sequences, columns=['ID', 'pA', 'pB'])


def merge_msas(file1, file2, file3):
    """
    Cette fonction lit trois fichiers MSA, les traite et les fusionne en un seul DataFrame.

    :param file1: Chemin vers le premier fichier FASTA
    :param file2: Chemin vers le deuxième fichier FASTA
    :param file3: Chemin vers le troisième fichier FASTA
    :return: DataFrame fusionné avec les colonnes ID, pA, et pB
    """
    # Lire chaque fichier MSA
    df1 = read_msa_fasta(file1)
    df2 = read_msa_fasta(file2)
    df3 = read_msa_fasta(file3)
    # Fusionner les DataFrames
    merged_df = pd.concat([df1, df2, df3], ignore_index=True)
    return merged_df



result_df = merge_msas('/Users/louiscarrel/Documents/Alignement_Project/largescale_kinase/DATA/ENSG00000107643/New_alignement/msa_s_exon_17_0.fasta', 
                       '/Users/louiscarrel/Documents/Alignement_Project/largescale_kinase/DATA/ENSG00000107643/New_alignement/msa_s_exon_17_1.fasta', 
                       '/Users/louiscarrel/Documents/Alignement_Project/largescale_kinase/DATA/ENSG00000107643/inter/undecided_sequences_17_0.fasta')


json_file_path = '/Users/louiscarrel/Documents/Alignement_Project/largescale_kinase/DATA/ENSG00000107643/Dict_species.json'
with open(json_file_path, 'r') as json_file:
    species_dict = json.load(json_file)

def add_species_column(df, species_mapping):
    """
    Ajoute une colonne 'Species' à un DataFrame en utilisant un dictionnaire de mappage.

    :param df: DataFrame original contenant les IDs des séquences
    :param species_mapping: Dictionnaire contenant le mappage ID -> Species
    :return: DataFrame modifié avec la colonne 'Species' ajoutée
    """
    df['Species'] = df['ID'].map(species_mapping).fillna('Unknown')
    return df

result_df = add_species_column(result_df, species_dict)
result_df.to_csv('/Users/louiscarrel/Documents/Alignement_Project/largescale_kinase/DATA/ENSG00000107643/result_with_species.csv', index=False)
result_df = result_df[result_df['Species'] != 'Unknown']
result_df.reset_index(drop=True, inplace=True)
result_df.to_csv('/Users/louiscarrel/Documents/Alignement_Project/largescale_kinase/DATA/ENSG00000107643/result_with_species_with_proba.csv', index=False)
print(result_df) 