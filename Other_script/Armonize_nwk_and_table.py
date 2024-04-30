import pandas as pd

# Chemins vers vos fichiers
path_csv = 'DATA/ENSG00000107643/filtered_species_for_pastml.tsv'
path_newick = 'DATA/ENSG00000107643/Species_list_time_tree.nwk'

# Fonction pour mettre en forme les noms d'espèces
def format_species_name(name):
    parts = name.lower().split()
    return '_'.join(part.capitalize() for part in parts)

# Lire le fichier CSV (ou TSV, ici on suppose séparé par des tabulations)
df = pd.read_csv(path_csv, sep='\t')

# Mettre à jour les noms d'espèces dans le DataFrame
df['Species'] = df['Species'].apply(format_species_name)

# Sauvegarder le fichier CSV modifié
df.to_csv('DATA/ENSG00000107643/filtered_species_for_pastml_formatted.tsv', index=False, sep='\t')

# Lire le fichier Newick
with open(path_newick, 'r') as file:
    newick_content = file.read()

# Mettre à jour les noms dans le fichier Newick
updated_newick = newick_content
for species in df['Species'].unique():
    original_species = ' '.join(part.capitalize() for part in species.split('_'))
    updated_newick = updated_newick.replace(original_species, species)

# Sauvegarder le fichier Newick modifié
with open('DATA/ENSG00000107643/Species_list_time_tree_formatted.nwk', 'w') as file:
    file.write(updated_newick)
