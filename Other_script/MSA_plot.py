# %%
import numpy as np
import pandas as pd
from pandas import json_normalize
from Bio import SeqIO
import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches




# %%
alph = ["a","c","d","e","f","g","h","i","k","l","m","n","p","q","r","s","t","v","w","y","-","x", "b", "z", "u"]
alph = [i.upper() for i in alph]

# %%
d_aa_id = dict(zip(alph, range(1,len(alph)+1)))

# %%
#convert seqMSA to num MSA
def toNum(msa, d_aa_id):
    Ls = len(msa[0])
    N = len(msa)
    #msaNum=np.zeros((N, Ls))
    msaNum = np.array([[d_aa_id[x]for x in seq] for seq in msa])
    return msaNum

# %%
def calculate_top_n_exons(msa, start_stop_df, num_exons=5):
    # ... (the rest of your function remains unchanged)

    # Calculate average lines for each exon and store in a list
    exon_lines = []
    for index, row in start_stop_df.iterrows():
        start_position = row['start']
        end_position = row['end']

        # Extract lines corresponding to the exon interval
        exon_lines.append(lines[:, start_position:end_position])

    # Calculate average sequence identity for each exon
    avg_sequence_identity = [np.nanmean(exon) for exon in exon_lines]

    # Sort exons based on average sequence identity
    sorted_exons = np.argsort(avg_sequence_identity)

    # Select the top 'num_exons' exons with the lowest average lines
    selected_exons = sorted_exons[:num_exons]

    return selected_exons, exon_lines

def plot_msa_v2(msa, d_aa_id, start_stop_df, path, sort_lines=True, dpi=100):
    seqQuery = msa[0]
    seqQueryNum =[d_aa_id[x] for x in seqQuery]
    #print(seqQuery, seqQueryNum)
    #seq = feature_dict["msa"][0]
    Ls = [len(seqQuery)]
    Ln = np.cumsum([0] + Ls) #cumulative length of sequences
    msaNum=toNum(msa,d_aa_id)
    gap = msaNum != 21             ## boulean indicating the presence of a gap ("-") in the MSA.
    qid = msaNum == seqQueryNum       ## boolean array where True indicates that the corresponding element in msa is equal to the query sequence (seq)
                                ## This comparison is done to identify positions in the MSA where the aligned sequences match the query sequence.
    #print(gap.shape, qid)
    gapid = np.stack([gap[:,Ln[i]:Ln[i+1]].max(-1) for i in range(len(Ls))],-1)
    '''This line creates a 2D array gapid that represents the presence of gaps within each segment defined by Ls. It checks for gaps in each segment and takes the maximum value across the columns.'''
    lines = []

    for g in np.unique(gapid, axis=0):
        i = np.where((gapid == g).all(axis=-1))
        qid_ = qid[i]
        gap_ = gap[i]
        seqid = np.stack([qid_[:,Ln[i]:Ln[i+1]].mean(-1) for i in range(len(Ls))],-1).sum(-1) / (g.sum(-1) + 1e-8)
        non_gaps = gap_.astype(float)
        non_gaps[non_gaps == 0] = np.nan
        if sort_lines:
            lines_ = non_gaps[seqid.argsort()] * seqid[seqid.argsort(),None]
        else:
            lines_ = non_gaps[::-1] * seqid[::-1,None]
        lines.append(lines_)
    
    
    lines = np.concatenate(lines,0)
    plt.figure(figsize=(8,5), dpi=dpi)
    plt.title("Sequence coverage")
    plt.imshow(lines,
              interpolation='nearest', aspect='auto',
              cmap="rainbow_r", vmin=0, vmax=1, origin='lower',
              extent=(0, lines.shape[1], 0, lines.shape[0]))
    for index, row in start_stop_df.iterrows():
        start_position = row['start']
        end_position = row['end']
        key_value = row['Key']

        # Draw vertical lines
        plt.axvline(x=start_position, color='gray', linestyle='--', linewidth=0.8)
        plt.axvline(x=end_position, color='gray', linestyle='--', linewidth=0.8)

        # Annotate the interval with the 'Key' value
        text_position = (start_position + end_position) / 2
        plt.text(text_position, 0, key_value, rotation=90, ha='center', va='bottom', color='black', fontsize=8)

   
    plt.plot((np.isnan(lines) == False).sum(0), color='black')
    plt.rcParams['savefig.facecolor']='white'
    plt.xlim(0,lines.shape[1])
    plt.ylim(0,lines.shape[0])
    plt.colorbar(label="Sequence identity to query")
    plt.xlabel("Positions")
    plt.ylabel("Sequences")
    plt.savefig(path, dpi=dpi, format='jpg', bbox_inches='tight')
    #plt.show()
    plt.close()

    coverage  = (np.isnan(lines) == False).sum(0)/len(msa)
    return coverage






def calculer_tangente(x, y):
    derivee = np.gradient(y, x)
    return derivee


def separate_start_end(df, value_column='Value'):
    """
    Sépare les colonnes 'start' et 'end' d'un dictionnaire dans une colonne spécifiée du DataFrame.

    """
    # Normalisation du dictionnaire dans la colonne spécifiée
    normalized_values = json_normalize(df[value_column].apply(eval))

    # Ajout des colonnes 'start' et 'end' au DataFrame existant
    df['start'] = normalized_values['start']
    df['end'] = normalized_values['end']

    # Suppression de la colonne d'origine si nécessaire
    df = df.drop(columns=[value_column])

    return df



if __name__ == "__main__":
# %%
    pp='P45983'
    dir_out = 'inter/'
    dir_GEMME = 'ENSG00000107643'
    pp_MSA = list(SeqIO.parse(dir_GEMME+f'/ali{pp}.fasta', "fasta"))
    l_seqMSA = [str(x.seq) for x in pp_MSA] 
    pathToSave = f'{dir_out}{pp}_MSA.jpg'
    start_stop = pd.read_csv('/Users/louiscarrel/Documents/Alignement_Project/inter/ENST00000374189_exon_coordinates_transcript.csv')
    start_stop_updated = separate_start_end(start_stop)

    lines= plot_msa_v2(l_seqMSA, d_aa_id, start_stop_updated, pathToSave, sort_lines=True, dpi=100)
    print(lines)


    

    x = np.arange(len(lines))
    tangente = calculer_tangente(x, lines)
    plt.plot(x, lines, label='Sequence identity to query')
    plt.plot(x, tangente, label='Tangente')
    plt.legend()
    plt.xlabel('Nucléotide')
    plt.ylabel('Coverage')
    plt.title('Tracé de la courbe et de la tangente')
    plt.show()
    