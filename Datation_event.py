import pandas as pd
from Bio import Phylo
import matplotlib.pyplot as plt

gene = "/Users/louiscarrel/Documents/Alignement_Project/largescale_kinase/DATA/ENSG00000107643/"
path_nwk = gene + 'result_with_species_with_proba.nwk'
path_tsv = gene + 'result_with_species_with_proba.tsv'
species_ref = "Homo_sapiens"


def get_distances(path_nwk, species_ref):
    # parse the tree file and builds the tree
    tree = Phylo.read(path_nwk, "newick")
    
    # get the names of the leaves
    leaves = {term.name for term in tree.get_terminals()}
    
    # Initialiser le dictionnaire pour stocker les distances
    d = {}
    
    # Calculer la distance pour chaque feuille par rapport à species_ref
    for l in leaves:
        # Empiriquement, il semble que nous devrions diviser par 2
        d[l.lower()] = tree.distance(species_ref, l) / 2
       
    return d





def plot_distances(path_nwk, path_tsv, species_ref):
    # Obtenir les distances
    distances = get_distances(path_nwk, species_ref)
    
    # Charger les données TSV
    data = pd.read_csv(path_tsv, sep='\t')  # Assurez-vous que le séparateur est correct
    data['Species'] = data['Species'].str.lower()  # Assurer la correspondance des cas
    
    # Ajouter les distances au DataFrame
    data['Distance'] = data['Species'].map(distances)
    
    # Supprimer les lignes où la distance est NaN (si aucune correspondance n'est trouvée)
    data.dropna(subset=['Distance'], inplace=True)
    
    # Tracer le graphique
    plt.figure(figsize=(10, 6))
    plt.scatter(data['Index'], data['Distance'], color='blue', alpha=0.5)
    plt.title('Distance from ' + species_ref + ' vs. probability to canonic path')
    plt.xlabel('probability')
    plt.ylabel('Distance')
    plt.grid(True)
    plt.show()

plot_distances(path_nwk, path_tsv, species_ref)
