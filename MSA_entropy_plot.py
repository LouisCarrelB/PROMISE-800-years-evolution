import numpy as np
from Bio import AlignIO
import matplotlib.pyplot as plt
import sys
import logomaker as lm
import pandas as pd
import Alignement_2_msa as ali
from Bio.Align import MultipleSeqAlignment
import os
import subprocess
from io import StringIO

def calculate_column_entropy(column):
    """Calculates Shannon entropy for a column of the alignment."""
    values, counts = np.unique(column, return_counts=True)
    probabilities = counts / counts.sum()
    entropy = -np.sum(probabilities * np.log2(probabilities, where=(probabilities!=0)))
    return entropy



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

def prepare_and_align(msa1, msa2):
    """Prepare and align two MSAs and return them separately."""
    concatenated_msa = msa1 + msa2  # Concatenate MSAs
    aligned_msa = align_msa_with_mafft(concatenated_msa)  # Align concatenated MSAs
    
    # Assuming both MSAs are of the same length after alignment
    msa1_aligned = MultipleSeqAlignment(aligned_msa[:len(msa1)])
    msa2_aligned = MultipleSeqAlignment(aligned_msa[len(msa1):])



    return msa1_aligned, msa2_aligned

def create_msa_logo(alignment):
    """
    Prend un objet MultipleSeqAlignment, extrait les séquences avec les gaps,
    et génère un DataFrame pour un logo plot.
    """
    import logomaker as lm

    # Garde les gaps dans les séquences
    seqs = [str(record.seq) for record in alignment]
    
    if len(seqs) == 0:
        print("Aucune séquence valide trouvée.")
        return None
    
    print(f'There are {len(seqs)} sequences, all of length {len(seqs[0])}')
    
    # Vérifie si toutes les séquences ont la même longueur
    if len(set(len(seq) for seq in seqs)) != 1:
        print("Erreur : Les longueurs de séquences ne sont pas uniformes.")
        return None

    # Convertit les séquences en DataFrame pour le logo
    msa_df = lm.alignment_to_matrix(sequences=seqs, to_type='counts', characters_to_ignore='X')
    
    return msa_df

def plot_msa_analysis(msa_path):
    """Loads an MSA, calculates entropy, generates a logo plot, and plots both vertically."""
    msa = AlignIO.read(msa_path, "fasta")
    entropies = calculate_msa_entropy(msa)
    msa_df = create_msa_logo(msa_path)
    
    if msa_df is None or msa_df.empty:
        return

    fig, axs = plt.subplots(2, 1, figsize=(12, 12), sharex=True)  # Change to 2 rows, 1 column with shared x-axis

    # Plot entropy histogram on top
    axs[0].bar(range(len(entropies)), entropies, color='blue')
    axs[0].set_ylabel('Entropy (bits)')
    axs[0].set_title(f"Entropy Histogram of MSA: {msa_path}")

    # Plot sequence logo on bottom
    lm.Logo(msa_df, ax=axs[1], fade_below=0., color_scheme='chemistry')
    axs[1].set_title(f"Sequence Logo of MSA: {msa_path}")
    axs[1].set_xlabel('Position in Alignment')

    plt.tight_layout()
    plt.show()


def plot_aa_distributions(msa_path, position):
    msa = list(AlignIO.read(msa_path, "fasta"))
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    aa_counts = {aa: 0 for aa in amino_acids}
    for record in msa:
        aa = record.seq[position]
        if aa in amino_acids:
            aa_counts[aa] += 1

    # Normalize the counts to get the probability distribution
    total = sum(aa_counts.values())
    aa_probs = {aa: count / total for aa, count in aa_counts.items()}

    # Plotting
    plt.bar(aa_probs.keys(), aa_probs.values())
    plt.title(f'Amino Acid Distribution at Position {position }')
    plt.xlabel('Amino Acid')
    plt.ylabel('Frequency')
    plt.show()


def plot_kl_distribution(msa1_path, msa2_path):
    """Charge les MSA, calcule la KL divergence pour chaque position, et affiche un histogramme."""
    msa1 = list(AlignIO.read(msa1_path, "fasta"))
    msa2 = list(AlignIO.read(msa2_path, "fasta"))
    
    kl_divergences = []
    
    # Assurez-vous que les deux MSA ont la même longueur
    if len(msa1[0].seq) != len(msa2[0].seq):
        raise ValueError("Les deux MSA doivent avoir la même longueur")
    
    # Calculer la divergence KL pour chaque position
    for i in range(len(msa1[0].seq)):
        col1 = [record.seq[i] for record in msa1 if record.seq[i] in 'ACDEFGHIKLMNPQRSTVWY']
        col2 = [record.seq[i] for record in msa2 if record.seq[i] in 'ACDEFGHIKLMNPQRSTVWY']
        
        dist1 = ali.calculate_aa_distribution(col1)
        dist2 = ali.calculate_aa_distribution(col2)
    
        
        kl1 = ali.calculate_kl_divergence(dist1, dist2)
        kl2 = ali.calculate_kl_divergence(dist2, dist1)
        
        kl_divergences.append((kl1 + kl2) / 2)  # Prenez la moyenne symétrique
        
    # Afficher l'histogramme
    plt.figure(figsize=(10, 6))
    plt.hist(kl_divergences, bins=30, color='blue', alpha=0.7)
    plt.title('Distribution de la divergence KL entre deux MSA')
    plt.xlabel('Divergence KL')
    plt.ylabel('Nombre de positions')
    plt.show()

def plot_combined_msa_analysis(msa1_path, msa2_path):
    """Loads two MSAs, calculates KL divergence, generates a logo plot for each, and plots KL divergence histogram."""
    # Load MSAs
    msa1 = list(AlignIO.read(msa1_path, "fasta"))
    msa2 = list(AlignIO.read(msa2_path, "fasta"))
    msa1,msa2 = prepare_and_align(msa1, msa2)



    # KL divergence calculations
    kl_divergences = []
    min_len = min(len(msa1[0].seq), len(msa2[0].seq))
    for i in range(min_len):
        col1 = [record.seq[i] for record in msa1 if record.seq[i] in 'ACDEFGHIKLMNPQRSTVWY']
        col2 = [record.seq[i] for record in msa2 if record.seq[i] in 'ACDEFGHIKLMNPQRSTVWY']
        
        dist1 = ali.calculate_aa_distribution(col1)
        dist2 = ali.calculate_aa_distribution(col2)
        
        kl1 = ali.calculate_kl_divergence(dist1, dist2)
        kl2 = ali.calculate_kl_divergence(dist2, dist1)
        kl_divergences.append((kl1 + kl2) / 2)
    
    # Create dataframes for logos
    msa_df2 = create_msa_logo(msa2)
    msa_df1 = create_msa_logo(msa1)
    

    # Create the figure with subplots
    fig = plt.figure(figsize=(18, 8))  # Width, Height
    gs = fig.add_gridspec(2, 2, width_ratios=[3, 1], hspace=0.4)  # 2 rows, 2 columns

    # Sequence logos
    ax_logo1 = fig.add_subplot(gs[0, 0])
    ax_logo2 = fig.add_subplot(gs[1, 0])
    if msa_df1 is not None:
        lm.Logo(msa_df1, ax=ax_logo1, fade_below=0., color_scheme='chemistry')
    if msa_df2 is not None:
        lm.Logo(msa_df2, ax=ax_logo2, fade_below=0., color_scheme='chemistry')

    # KL divergence histogram
    ax_kl = fig.add_subplot(gs[:, 1])  # Span all rows in the second column
    ax_kl.bar(range(len(kl_divergences)), kl_divergences, color='blue')
    ax_kl.set_ylabel('KL Divergence (bits)')
    ax_kl.set_title('KL Divergence Histogram between MSAs')
    plt.show()


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script_name.py path_to_msa1_file.fasta, path_to_msa2_file.fasta")
    else:
        msa1_path = sys.argv[1]
        msa2_path = sys.argv[2]
        plot_msa_analysis(msa1_path)
        plot_combined_msa_analysis(msa1_path,msa2_path)
