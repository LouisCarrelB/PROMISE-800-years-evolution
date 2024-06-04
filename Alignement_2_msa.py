from Bio import AlignIO
import numpy as np

from Bio import SeqIO, AlignIO
import numpy as np
import math
import tqdm
import matplotlib.pyplot as plt
import requests
import pandas as pd

from Bio.Align import MultipleSeqAlignment

def calculate_presence_matrix(alignment, special_positions, multiplier=5):
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'  # Liste des 20 acides aminés standards
    num_positions = len(alignment[0])  # Longueur de la première séquence, suppose que toutes ont la même longueur
    presence_matrix = np.zeros((len(amino_acids), num_positions))

    # Calculer le ratio de présence pour chaque acide aminé à chaque position
    for pos in range(num_positions):
        column = [record.seq[pos] for record in alignment if record.seq[pos] in amino_acids]
        total_aas = len(column)

        if total_aas > 0:  # Éviter la division par zéro
            aa_counts = {aa: column.count(aa) for aa in amino_acids}
            for i, aa in enumerate(amino_acids):
                presence_matrix[i, pos] = aa_counts[aa] / total_aas
                # Appliquer le multiplicateur si la position est spéciale
                if pos in special_positions:
                    presence_matrix[i, pos] *= multiplier

    return presence_matrix

def plot_presence_matrix(presence_matrix, amino_acids):
    fig, ax = plt.subplots(figsize=(20, 10))  # Vous pouvez ajuster la taille selon vos besoins
    cax = ax.matshow(presence_matrix, cmap='viridis')  # Choix de la colormap
    plt.title("Acide Aminé Presence Matrix")
    plt.xlabel("Position in Alignment")
    plt.ylabel("Acide Aminé")
    
    # Ajouter une colorbar
    fig.colorbar(cax)

    # Définir les ticks de l'axe y pour montrer les acides aminés
    ax.set_yticks(np.arange(len(amino_acids)))
    ax.set_yticklabels(amino_acids)

    # Afficher les ticks de l'axe x pour chaque position avec un intervalle personnalisé
    ax.set_xticks(np.arange(0, presence_matrix.shape[1], 10))
    ax.set_xticklabels(np.arange(1, presence_matrix.shape[1] + 1, 10))

    plt.show()


def calculate_sequence_score(sequence, presence_matrix, amino_acids='ACDEFGHIKLMNPQRSTVWY'):
    score = 0
    for pos, aa in enumerate(sequence):
        if aa in amino_acids:
            aa_index = amino_acids.index(aa)
            score += presence_matrix[aa_index, pos]
    return score

def calculate_entropy(acid_set, total_count):
    entropy = 0
    for acid in acid_set:
        p = acid_set[acid] / total_count
        entropy -= p * np.log2(p)
    return entropy


def find_mutually_exclusive_positions(msa1, msa2):
    if len(msa1[0].seq) != len(msa2[0].seq):
        raise ValueError("Les deux MSA doivent avoir la même longueur d'alignement")
    
    mutually_exclusive_positions = []
    alignment_length = len(msa1[0].seq)

    for i in range(alignment_length):
        count_set1 = {seq.seq[i]: 0 for seq in msa1}
        count_set2 = {seq.seq[i]: 0 for seq in msa2}

        for seq in msa1:
            count_set1[seq.seq[i]] += 1
        for seq in msa2:
            count_set2[seq.seq[i]] += 1

        set1 = set(count_set1.keys())
        set2 = set(count_set2.keys())

        total1 = sum(count_set1.values())
        total2 = sum(count_set2.values())

        entropy1 = calculate_entropy(count_set1, total1)
        entropy2 = calculate_entropy(count_set2, total2)

        # Vérifier si les ensembles sont disjoints et prendre en compte l'entropie
        if set1.isdisjoint(set2) and (entropy1 > 0 or entropy2 > 0):
            mutually_exclusive_positions.append((i, entropy1, entropy2))

    return mutually_exclusive_positions



def align_sequences_to_msas(input_file, msa1_path, msa2_path, output_dir):
    # Charger les MSA
    msa1 = AlignIO.read(msa1_path, "fasta")
    msa2 = AlignIO.read(msa2_path, "fasta")
    sequences = list(SeqIO.parse(input_file, "fasta"))

    # Calculer les matrices de présence pour chaque MSA
    special_positions = []  # Définir ou trouver ces positions
    presence_matrix_msa1 = calculate_presence_matrix(msa1, special_positions, 5)
    presence_matrix_msa2 = calculate_presence_matrix(msa2, special_positions, 5)

    can, alt, undecided = [], [], []

    for sequence in sequences:
        score_msa1 = calculate_sequence_score(sequence.seq, presence_matrix_msa1)
        score_msa2 = calculate_sequence_score(sequence.seq, presence_matrix_msa2)
        
        pA = math.exp(score_msa1) / (math.exp(score_msa1) + math.exp(score_msa2))
        pB = math.exp(score_msa2) / (math.exp(score_msa1) + math.exp(score_msa2))

        # Créer un nouveau SeqRecord avec les scores dans la description
        new_seq_record = sequence
        new_seq_record.description = f"pA: {pA:.3f}, pB: {pB:.3f}"

        if pA > 0.7:  # Ajouter la séquence à msa1 si pA est significativement plus grand
            can.append(new_seq_record)
        elif pB > 0.7:  # Ajouter la séquence à msa2 si pB est significativement plus grand
            alt.append(new_seq_record)
        else:
            undecided.append(new_seq_record)  # Séquences indécises

    # Créer des objets MultipleSeqAlignment
    can_alignment = MultipleSeqAlignment(can)
    alt_alignment = MultipleSeqAlignment(alt)

    # Sauvegarder les nouveaux MSA et les séquences indécises
    AlignIO.write(can_alignment, f"{output_dir}/can_ali.fasta", "fasta")
    AlignIO.write(alt_alignment, f"{output_dir}/alt_ali.fasta", "fasta")
    SeqIO.write(undecided, f"{output_dir}/undecided_sequences.fasta", "fasta")
    




def fetch_species_name(gene_id):
    # Interroger l'API d'Ensembl pour obtenir le nom de l'espèce à partir de l'ID de gène
    server = "https://rest.ensembl.org"
    ext = f"/lookup/id/{gene_id}?content-type=application/json"
    response = requests.get(server+ext)
    if response.status_code == 200:
        data = response.json()
        return data.get('species', 'Unknown')
    else:
        return 'Unknown'
    
def create_gene_species_dict(exon_table_path):
    # Charger la table d'exons
    exon_df = pd.read_csv(exon_table_path)
    # Créer un dictionnaire sans doublons où GeneID mappe à Species
    gene_species_dict = pd.Series(exon_df.Species.values, index=exon_df.GeneID).to_dict()
    return gene_species_dict

def classify_genes(msa1_path, msa2_path, msa3_path, exon_table_path):
    # Créer le dictionnaire GeneID à Species
    gene_species_dict = create_gene_species_dict(exon_table_path)

    # Charger les séquences des fichiers FASTA
    can_genes = {seq.id: 'can' for seq in SeqIO.parse(msa1_path, "fasta")}
    alt_genes = {seq.id: 'alt' for seq in SeqIO.parse(msa2_path, "fasta")}
    undecided_genes = {seq.id: 'undecided' for seq in SeqIO.parse(msa3_path, "fasta")}

    # Dictionnaire pour maintenir les classifications des gènes
    gene_classification = {}

    # Classer chaque gène et associer le nom de l'espèce
    unique_gene_ids = set(can_genes).union(alt_genes, undecided_genes)
    for gene_id in unique_gene_ids:
        classes = []
        if gene_id in can_genes:
            classes.append(can_genes[gene_id])
        if gene_id in alt_genes:
            classes.append(alt_genes[gene_id])
        if gene_id in undecided_genes:
            classes.append(undecided_genes[gene_id])
        
        # Récupérer le nom de l'espèce depuis le dictionnaire
        species_name = gene_species_dict.get(gene_id, 'Unknown')

        # Stocker les informations dans le dictionnaire
        gene_classification[gene_id] = {
            'Species': species_name,
            'Group': '/'.join(classes)
        }

    # Créer le DataFrame
    gene_df = pd.DataFrame.from_dict(gene_classification, orient='index')
    return gene_df



if __name__ == "__main__":
    input_file = "ENSG00000107643/thoraxe_2/non_kept_sequences.fasta"
    msa1_path = "ENSG00000107643/thoraxe_2/can.fasta"
    msa2_path = "ENSG00000107643/thoraxe_2/alt.fasta"
    output_dir = "ENSG00000107643/thoraxe_2/aligned_sequences"
    align_sequences_to_msas(input_file, msa1_path, msa2_path, output_dir)


    msa1_path = f"{output_dir}/can_ali.fasta"
    msa2_path = f"{output_dir}/alt_ali.fasta"
    undecided_path = f"{output_dir}/undecided_sequences.fasta"

    exon_table_path = "ENSG00000107643/thoraxe_2/s_exon_table.csv"
    

    gene_df = classify_genes(msa1_path, msa2_path, undecided_path,exon_table_path)
    print(gene_df)
    gene_df.to_csv('gene_group_classification.csv', index=True)
    species_only_df = gene_df[['Species']]
    species_only_df.to_csv('species_only.csv', index=False)

    # Sauvegarde avec les colonnes 'Species' et 'Group', sans les noms de ligne
    gene_df['Species'] = gene_df['Species'].str.lower().str.capitalize().str.replace(' ', '_')
    species_to_remove = ['Oryzias_melastigma', 'Salmo_trutta']
    gene_classification_table = gene_df[~gene_df['Species'].isin(species_to_remove)]
    species_group_df = gene_df[['Species', 'Group']]


    species_group_df.to_csv('species_and_group.tsv', sep='\t', index=False)