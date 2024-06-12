import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist
import argparse
from tqdm import tqdm 

def lire_a3m(fichier_a3m):
    with open(fichier_a3m, 'r') as file:
        sequences = [line.strip() for line in file if not line.startswith('>')]
    return sequences

def trouver_insertions(sequences):
    insertions = []
    for seq in sequences:
        parts = seq.split('-')
        insertions.extend([part for part in parts if len(part) > 0])
    return insertions


def cluster_insertions(insertions, seuil=0.2):
    n = len(insertions)
    # Initialiser une matrice de distances à zéro
    distances = np.zeros((n * (n - 1)) // 2)
    k = 0
    for i in tqdm(range(n), desc="Calculating distances"):
        for j in range(i + 1, n):
            # Calculer la distance comme le rapport de différences sur la longueur de la plus longue insertion
            diff_count = sum(1 for a, b in zip(insertions[i], insertions[j]) if a != b)
            diff_count += abs(len(insertions[i]) - len(insertions[j])) # Compter les caractères supplémentaires comme différences
            max_len = max(len(insertions[i]), len(insertions[j]))
            distances[k] = diff_count / max_len
            k += 1

    Z = linkage(distances, 'ward')
    clusters = fcluster(Z, seuil, criterion='distance')
    return clusters, insertions


def filtrer_par_taille(clusters, insertions, taille_min=5):
    cluster_dict = {}
    for cluster_id, insertion in tqdm(zip(clusters, insertions),desc="Filtering by size"):
        if len(insertion) >= taille_min:
            if cluster_id in cluster_dict:
                cluster_dict[cluster_id].append(insertion)
            else:
                cluster_dict[cluster_id] = [insertion]
    return cluster_dict

def main():
    parser = argparse.ArgumentParser(description="Trouver des insertions courantes dans un fichier .a3m.")
    parser.add_argument("chemin_fichier_a3m", type=str, help="Path to .a3m file")
    args = parser.parse_args()

    sequences = lire_a3m(args.chemin_fichier_a3m)
    insertions = trouver_insertions(sequences)
    clusters, insertions = cluster_insertions(insertions)
    clusters_filtrés = filtrer_par_taille(clusters, insertions)                                                                     
    for cluster_id, ins in clusters_filtrés.items():
        print(f"Cluster {cluster_id}: {ins}")

if __name__ == "__main__":
    main()
