import editdistance
import re 
import glob
from Bio import SeqIO, AlignIO
import requests
import csv 
import pandas as pd 
import sys 

def read_pir_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    return lines

def find_sequence_in_pir(sequence, pir_lines, max_distance=20):
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
                # Extraire la partie entre le premier et le prochain espace
                space_index = current_sequence.find(' ')
                match = re.search(r';([^ ]+)\s([^ ]+)', current_sequence)
                if match:
                    first_identifier.append(match.group(1))
                    second_identifier.append(match.group(2))
    return first_identifier,second_identifier,found_sequences


def extract_uni_id(sequence_id):
    # Expression régulière pour extraire ce qui se trouve après le "_"
    pattern = re.compile(r"_(.+)")
    # Recherche de la correspondance dans sequence_id
    match = pattern.search(sequence_id)
    # Si une correspondance est trouvée, retourner le résultat après le "_"
    if match:
        return match.group(1)
    else:
        # Si aucun "_" n'est trouvé, retourner l'entier de sequence_id
        try:
            return sequence_id
        except ValueError:
            # Si la conversion en entier échoue, retourner None
            return None
        
def get_organism_name(uniprot_id):
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
    # Charger le fichier CSV en spécifiant le délimiteur
    s_exons_df = pd.read_csv(file_path, delimiter=',')

    # Assurez-vous que le nom de la colonne 'GeneID' est correctement ciblé
    unique_gene_ids = s_exons_df['GeneID'].astype(str).unique()

    # Créer un dictionnaire pour lier une Species à GeneID
    gene_species_dict = {}
    for gene_id in unique_gene_ids:
        gene_species_dict[gene_id] = list(s_exons_df[s_exons_df['GeneID'] == gene_id]['Species'].unique())

    return gene_species_dict



if __name__ == "__main__":
    csv_data = []
    if len(sys.argv) < 2:
        print("Usage: python Alignement_PIR.py <gene_name>")
        sys.exit(1)
        
    gene_name = sys.argv[1]


    GENE = "DATA/" + gene_name + "/"


    
    pir_file_path = GENE + 'thoraxe/phylosofs/transcripts.pir'  
    s_exons = GENE + "thoraxe/s_exon_table.csv"


    for a3m_fichier in glob.glob(GENE + "*.a3m"):
            sequences_a3m = list(SeqIO.parse(a3m_fichier, "fasta"))[:30]
            for sequence_a3m in sequences_a3m:
                pattern = re.compile(r"_(.+)")

                target_sequence = sequence_a3m.seq
                sequence_id = sequence_a3m.id

                print(extract_uni_id(sequence_id))
                organism_name = get_organism_name(extract_uni_id(sequence_id))
                pir_lines = read_pir_file(pir_file_path)
                results_ID_g,results_ID_t,results_tot = find_sequence_in_pir(target_sequence, pir_lines)
                print("pour la séquence : ",sequence_id)
                print ("organisme : " ,get_organism_name(extract_uni_id(sequence_id)))
                if results_ID_g:
                    
                    print("Séquences trouvées :")
                    for sequence, distance in results_tot:
                        print(f"{sequence} (Distance d'édition : {distance})")
                    print("Identifiants du transcrit trouvés :")
                    for identifier in results_ID_t:
                        print(identifier)
                    print("Identifiants du gène trouvés :")
                    for identifier in results_ID_g:
                        print(identifier)
                    for i, (sequence, distance) in enumerate(results_tot):
                        print(f"{sequence} (Distance d'édition : {distance}")
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
                    print("Aucun identifiant correspondant trouvé.") 






    dic_ensembl = create_gene_species_dict(s_exons)


    csv_data_df = pd.DataFrame(csv_data)

    # Ajouter une nouvelle colonne 'organism_pir' en utilisant le dictionnaire
    csv_data_df['organism_pir'] = csv_data_df['Gene IDs'].map(dic_ensembl)

    # Afficher le DataFrame mis à jour
    print(csv_data_df)

    csv_data_df.to_csv(GENE +'inter/a3m_to_PIR.csv', index=False)
        
    

    