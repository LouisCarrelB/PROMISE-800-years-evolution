from Bio import pairwise2
from Bio.SubsMat import MatrixInfo
from sklearn.preprocessing import OneHotEncoder
from sklearn.cluster import KMeans
from sklearn.preprocessing import normalize
import numpy as np 

# Chargement de la matrice BLOSUM62
blosum62 = MatrixInfo.blosum62
print(blosum62)
def aligner_sequence(sequence, msa):
    best_score = float("-inf")
    best_alignment = None

    for seq in msa:
        alignments = pairwise2.align.localds(seq, sequence, blosum62, -10, -0.5)
        for alignment in alignments:
            score = alignment[2]
            if score > best_score:
                best_score = score
                best_alignment = alignment

    return best_alignment

def recalculer_profil_hmm(msa):
    max_len = max(len(seq) for seq in msa)
    msa_aligned = [seq.ljust(max_len, '-') for seq in msa]
    
    # Conversion des séquences en vecteurs one-hot
    encoder = OneHotEncoder(dtype=int)
    msa_enc = encoder.fit_transform([[c for c in seq] for seq in msa_aligned])
    
    # Clustering des vecteurs one-hot pour obtenir les émissions probables
    kmeans = KMeans(n_clusters=2, random_state=0)
    kmeans.fit(msa_enc.toarray())  # Convertir la matrice creuse en tableau pour KMeans
    emissions = normalize(kmeans.cluster_centers_, norm='l1')
    
    return emissions







'''

Emissions probables par état : Les émissions probables par état sont des matrices qui indiquent
 les probabilités de chaque symbole (ou caractère) apparaissant dans chaque position d'une séquence, pour chaque état du modèle HMM. Dans notre exemple, 
 chaque ligne de la matrice `emissions` correspond à un état, et chaque colonne correspond à un symbole possible. Les valeurs dans la matrice représentent les probabilités que chaque symbole soit émis par chaque état. 
 Une valeur élevée pour un symbole dans un état indique que ce symbole a une forte probabilité d'être émis par cet état.



'''


msa_1 = [
    "ACGT",
    "ACTT",
    "ACGT",
    "ACCT"
]

msa_2 = [
    "TACG",
    "TTCC",
    "TGTT",
    "TAGT"
]

nouvelle_sequence = "ACGGT"

best_alignment_msa_1 = aligner_sequence(nouvelle_sequence, msa_1)
best_alignment_msa_2 = aligner_sequence(nouvelle_sequence, msa_2)
msa_1.append(nouvelle_sequence)
msa_2.append(nouvelle_sequence)
profil_hmm_msa_1 = recalculer_profil_hmm(msa_1)
profil_hmm_msa_2 = recalculer_profil_hmm(msa_2)

print("Profil HMM pour le MSA 1 mis à jour :")
print("Emissions probables par état :", profil_hmm_msa_1)

print("Profil HMM pour le MSA 2 mis à jour :")
print("Emissions probables par état :", profil_hmm_msa_2)


# Calcul de la cohérence des profils HMM pour chaque MSA
def evaluer_cohérence_profil(profil_hmm):
    # Ici, nous calculons la cohérence en sommant les éléments des matrices d'émissions et de transition
    coherence = np.sum(profil_hmm[0]) + np.sum(profil_hmm[1])
    return coherence

# Évaluer la cohérence des profils HMM pour chaque MSA
coherence_msa_1 = evaluer_cohérence_profil(profil_hmm_msa_1)
coherence_msa_2 = evaluer_cohérence_profil(profil_hmm_msa_2)

# Comparaison des résultats de cohérence
print("Cohérence du MSA 1 :", coherence_msa_1)
print("Cohérence du MSA 2 :", coherence_msa_2)

# Choix du MSA le plus cohérent
if coherence_msa_1 > coherence_msa_2:
    print("Le MSA 1 est plus cohérent.")
    msa_coherent = msa_1
    profil_hmm_coherent = profil_hmm_msa_1
else:
    print("Le MSA 2 est plus cohérent.")
    msa_coherent = msa_2
    profil_hmm_coherent = profil_hmm_msa_2

# Vous pouvez ensuite utiliser msa_coherent et profil_hmm_coherent pour votre analyse ultérieure
