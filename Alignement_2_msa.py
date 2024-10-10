from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd
from scipy.special import rel_entr
from Bio.Align import MultipleSeqAlignment
from Bio.Align import MultipleSeqAlignment
import os
import subprocess
from io import StringIO
import sys 


def align_msa_with_mafft(msa_records):
    """Align MSA using MAFFT."""
    # Convert BioPython records into a FASTA formatted string
    fasta_format = "".join(f">{record.id}\n{str(record.seq)}\n" for record in msa_records)
    
    # Set up the MAFFT command
    mafft_cmd = ["mafft", "--auto", "--quiet", "-"]
    
    # Execute MAFFT using subprocess with stdin and stdout
    process = subprocess.Popen(mafft_cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = process.communicate(fasta_format)
    
    if process.returncode != 0:
        raise Exception(f"MAFFT error: {stderr}")
    
    # Read the aligned sequences from MAFFT output
    return AlignIO.read(StringIO(stdout), "fasta")

def prepare_and_align_three_msas(msa1, msa2, msa3):
    """Prepare and align three MSAs and return them separately."""
    concatenated_msa = msa1 + msa2 + msa3  # Concatenate MSAs
    aligned_msa = align_msa_with_mafft(concatenated_msa)  # Align concatenated MSAs
    
    # Assuming all MSAs are of the same length after alignment
    msa1_aligned = MultipleSeqAlignment(aligned_msa[:len(msa1)])
    msa2_aligned = MultipleSeqAlignment(aligned_msa[len(msa1):len(msa1)+len(msa2)])
    msa3_aligned = MultipleSeqAlignment(aligned_msa[len(msa1)+len(msa2):])
    
    return msa1_aligned, msa2_aligned, msa3_aligned

def calculate_presence_matrix(alignment, position_weights):
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'  # Liste des 20 acides aminés standards
    num_positions = len(alignment[0].seq)  # Longueur de la première séquence, suppose que toutes ont la même longueur
    presence_matrix = np.zeros((len(amino_acids), num_positions))

    # Calculer le ratio de présence pour chaque acide aminé à chaque position
    for pos in range(num_positions):
        column = [record.seq[pos] for record in alignment if record.seq[pos] in amino_acids]
        total_aas = len(column)

        if total_aas > 0:  # Éviter la division par zéro
            aa_counts = {aa: column.count(aa) for aa in amino_acids}
            for i, aa in enumerate(amino_acids):
                presence_matrix[i, pos] = aa_counts[aa] / total_aas
                # Appliquer un poids spécifique si la position est jugée importante
                if pos in position_weights:
                    presence_matrix[i, pos] *=  position_weights[pos]
                else:
                    presence_matrix[i, pos] *= 1  # Si aucune information de poids, multiplier par 1 (ne change rien)

    return presence_matrix

def plot_presence_matrix(presence_matrix, amino_acids='ACDEFGHIKLMNPQRSTVWY'):
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



def align_sequences_to_msas(input_file, msa1_path, msa2_path, output_dir, iterations=1):
    # Charger les MSA
    msa1_ori = list(AlignIO.read(msa1_path, "fasta"))
    msa2_ori = list(AlignIO.read(msa2_path, "fasta"))
    sequences = list(SeqIO.parse(input_file, "fasta"))

    msa1_ori, msa2_ori, sequences = prepare_and_align_three_msas(msa1_ori, msa2_ori, sequences)
    msa1_ori = list(msa1_ori)
    msa2_ori = list(msa2_ori)
    sequences = list(sequences)
    for i in range(iterations):
        if i == 0:
            # Calculer les positions importantes spécifiques à chaque MSA à chaque itération
            special_positions = find_important_positions_with_weights(msa1_ori, msa2_ori)
            # Calculer les matrices de présence pour chaque MSA
            presence_matrix_msa1 = calculate_presence_matrix(msa1_ori, special_positions)
            presence_matrix_msa2 = calculate_presence_matrix(msa2_ori, special_positions)
            new_msa1 = []
            new_msa2 = []
            new_msa1 = msa1_ori.copy()
            new_msa2 = msa2_ori.copy()
        else :
            # Calculer les positions importantes spécifiques à chaque MSA à chaque itération
            special_positions = find_important_positions_with_weights(msa1, msa2)
            # Calculer les matrices de présence pour chaque MSA
            presence_matrix_msa1 = calculate_presence_matrix(msa1, special_positions)
            presence_matrix_msa2 = calculate_presence_matrix(msa2, special_positions)
            new_msa1 = []
            new_msa2 = []
            new_msa1 = msa1_ori.copy()
            new_msa2 = msa2_ori.copy()
       


        undecided = []
        for sequence in sequences:
            print(sequence)
            score_msa1 = calculate_sequence_score(sequence.seq, presence_matrix_msa1)
            score_msa2 = calculate_sequence_score(sequence.seq, presence_matrix_msa2)

            pA = math.exp(score_msa1) / (math.exp(score_msa1) + math.exp(score_msa2))
            pB = math.exp(score_msa2) / (math.exp(score_msa1) + math.exp(score_msa2))

            new_seq_record = sequence
            new_seq_record.description = f"pA: {pA:.3f}, pB: {pB:.3f}"

            if pA > 0.7:
                new_msa1.append(new_seq_record)
            elif pB > 0.7:
                new_msa2.append(new_seq_record)
            else:
                undecided.append(new_seq_record)

        # Update MSAs
        msa1 = new_msa1
        msa2 = new_msa2

    # Sauvegarder les nouveaux MSA et les séquences indécises
    AlignIO.write(MultipleSeqAlignment(msa1), f"{output_dir}/can_ali.fasta", "fasta")
    AlignIO.write(MultipleSeqAlignment(msa2), f"{output_dir}/alt_ali.fasta", "fasta")
    SeqIO.write(undecided, f"{output_dir}/undecided_sequences.fasta", "fasta")


def create_gene_species_dict(exon_table_path):
    # Charger la table d'exons
    exon_df = pd.read_csv(exon_table_path)
    # Créer un dictionnaire sans doublons où GeneID mappe à Species
    gene_species_dict = pd.Series(exon_df.Species.values, index=exon_df.GeneID).to_dict()
    return gene_species_dict

def classify_genes(msa1, msa2, msa3, exon_table_path):
    # Créer le dictionnaire GeneID à Species
    gene_species_dict = create_gene_species_dict(exon_table_path)

    # Créer les séquences des gènes
    can_genes = {seq.id: 'can' for seq in msa1}
    alt_genes = {seq.id: 'alt' for seq in msa2}
    undecided_genes = {seq.id: 'undecided' for seq in msa3}

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

def calculate_kl_divergence(p, q):
    """Calcule la divergence KL de p vers q en utilisant rel_entr pour gérer les probabilités de zéro."""
    p = np.asarray(p, dtype=np.float64)
    q = np.asarray(q, dtype=np.float64)
    # Ajout d'une petite constante pour éviter le log de zéro
    epsilon = 1e-10
    p += epsilon
    q += epsilon
    return np.sum(rel_entr(p, q))

def calculate_aa_distribution(column):
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    counts = {aa: 0 for aa in amino_acids}
    for aa in column:
        if aa in amino_acids:
            counts[aa] += 1
    total = sum(counts.values())
    distribution = [counts[aa] / total if total > 0 else 0 for aa in amino_acids]
    return distribution

def find_important_positions_with_weights(msa1, msa2):
    important_positions = {}
    all_kl1 = []
    all_kl2 = []
    for i in range(len(msa1[0].seq)):
        print(record.seq[i] for record in msa1)
        col1 = [record.seq[i] for record in msa1]
        col2 = [record.seq[i] for record in msa2]
        
        dist1 = calculate_aa_distribution(col1)
        dist2 = calculate_aa_distribution(col2)
        
        kl1 = calculate_kl_divergence(dist1, dist2)
        kl2 = calculate_kl_divergence(dist2, dist1)
        
        all_kl1.append(kl1)
        all_kl2.append(kl2)
        
        # Utilisation de la moyenne normalisée des divergences de KL comme poids
        mean_kl = (kl1 + kl2) / 2
        if mean_kl > 0:
            important_positions[i] = mean_kl
    
  

    return important_positions



if __name__ == "__main__":
    input_file = sys.argv[1]

    can = sys.argv[2]
    alt = sys.argv[3]





    output_dir = sys.argv[4]

    msa1 = list(AlignIO.read(can, "fasta"))
    msa2 = list(AlignIO.read(alt, "fasta"))


    
    align_sequences_to_msas(input_file, can, alt, output_dir)
