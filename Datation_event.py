import pandas as pd
from Bio import Phylo
import matplotlib.pyplot as plt
import statsmodels.api as sm
import numpy as np

gene = "/Users/louiscarrel/Documents/Alignement_Project/largescale_kinase/DATA/ENSG00000107643/"
path_nwk = gene + 'Species_list_time_tree.nwk'
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


def plot_probability_curve_with_annotations(path_nwk, path_tsv, species_ref):
    """
    Plot a probability curve with annotations, grouping data by evolutionary distance and averaging probabilities.
    
    Args:
    path_nwk (str): Path to the Newick format tree file.
    path_tsv (str): Path to the TSV file containing species and probability data.
    species_ref (str): Reference species name used for distance calculation.
    """
    distances = get_distances(path_nwk, species_ref)
    data = pd.read_csv(path_tsv, sep='\t')
    data['Species'] = data['Species'].str.lower()
    data['Distance'] = data['Species'].map(distances)
    data.dropna(subset=['Distance'], inplace=True)
    
    plt.figure(figsize=(12, 8))
    plt.scatter(data['Distance'], data['Index'], alpha=0.5, label='Data Points')
    
    # Add a lowess smoothed curve
    lowess = sm.nonparametric.lowess
    smoothed = lowess(data['Index'], data['Distance'], frac=0.3)
    plt.plot(smoothed[:, 0], smoothed[:, 1], color='red', label='Lowess Smoothing')
    
    # Annotate every 10th point for clarity
    for index, row in data.iterrows():
        if index % 60 == 0:
            plt.annotate(row['Species'], (row['Distance'], row['Index']), textcoords="offset points", xytext=(0,10), ha='center')
    
    plt.title('Probability vs. Evolutionary Distance with Annotations for the canonique exon')
    plt.xlabel('Evolutionary Distance')
    plt.ylabel('Probability')
    plt.grid(True)
    plt.legend()
    plt.show()




def plot_grouped_histogram_by_distance(path_nwk, path_tsv, species_ref):
    # Get distances
    distances = get_distances(path_nwk, species_ref)
    
    # Load TSV data
    data = pd.read_csv(path_tsv, sep='\t')
    data['Species'] = data['Species'].str.lower()
    data['Distance'] = data['Species'].map(distances)
    data.dropna(subset=['Distance'], inplace=True)
    
    # Define probability bins and labels
    probability_bins = [0, 0.25, 0.5, 0.75, 1.0]
    probability_labels = ['0-25%', '25-50%', '50-75%', '75-100%']
    data['Probability Group'] = pd.cut(data['Index'], bins=probability_bins, labels=probability_labels, include_lowest=True)
    
    # Define the distance range groups. You may want to adjust these based on your dataset.
    distance_bins = pd.qcut(data['Distance'], q=250, precision=2, duplicates='drop')  # Quantile-based to ensure even distribution
    data['Distance Group'] = distance_bins
    
    # Plotting
    grouped_data = data.groupby(['Distance Group', 'Probability Group']).size().unstack(fill_value=0)
    grouped_data.plot(kind='bar', stacked=False, figsize=(12, 8))
    
    plt.title('Probability Distribution Across Distance Groups for canonic path')
    plt.xlabel('Distance Groups')
    plt.ylabel('Count of Species')
    plt.legend(title='Probability Percentages')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()





def plot_distance_histogram_by_probability_groups(path_nwk, path_tsv, species_ref):
    # Get distances
    distances = get_distances(path_nwk, species_ref)
    
    # Load TSV data
    data = pd.read_csv(path_tsv, sep='\t')
    data['Species'] = data['Species'].str.lower()
    data['Distance'] = data['Species'].map(distances)
    data.dropna(subset=['Distance'], inplace=True)
    
    # Define finer probability bins and labels for every 5%
    probability_bins = list(np.arange(0, 1.05, 0.05))
    probability_labels = [f'{int(100 * i)}-{int(100 * (i + 0.05))}%' for i in np.arange(0, 1, 0.05)]
    data['Probability Group'] = pd.cut(data['Index'], bins=probability_bins, labels=probability_labels, include_lowest=True)
    
    # Define the distance range groups
    distance_bins = pd.qcut(data['Distance'], q=500, precision=2, duplicates='drop')
    data['Distance Group'] = distance_bins
    
    # Obtain a color map
    cmap = plt.get_cmap('viridis')  # You can change 'viridis' to any other colormap, like 'plasma' or 'rainbow'
    colors = cmap(np.linspace(0, 1, len(probability_labels)))

    # Plotting
    grouped_data = data.groupby(['Distance Group', 'Probability Group']).size().unstack(fill_value=0)
    grouped_data.plot(kind='bar', stacked=True, figsize=(12, 8), color=colors)
    
    plt.title('Probability Percentage Groups Across Distance Groups')
    plt.xlabel('Distance Groups')
    plt.ylabel('Count of Species')
    plt.legend(title='Probability Percentages', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks(rotation=60)
    plt.tight_layout()
    plt.show()



def plot_probability_histogram_by_distance_groups(path_nwk, path_tsv, species_ref):
    # Get distances
    distances = get_distances(path_nwk, species_ref)
    
    # Load TSV data
    data = pd.read_csv(path_tsv, sep='\t')
    data['Species'] = data['Species'].str.lower()
    data['Distance'] = data['Species'].map(distances)
    data.dropna(subset=['Distance'], inplace=True)
    
    # Define finer probability bins and labels for every 5%
    probability_bins = list(np.arange(0, 1.05, 0.05))
    probability_labels = [f'{int(100 * i)}-{int(100 * (i + 0.05))}%' for i in np.arange(0, 1, 0.05)]
    data['Probability Group'] = pd.cut(data['Index'], bins=probability_bins, labels=probability_labels, include_lowest=True)
    
    # Define the distance range groups and extract categories for legend
    distance_bins = pd.qcut(data['Distance'], q=50, precision=2, duplicates='drop')
    data['Distance Group'] = distance_bins
    distance_categories = distance_bins.cat.categories  # Extract categories for legend
    
    # Obtain a color map
    cmap = plt.get_cmap('viridis')
    colors = cmap(np.linspace(0, 1, len(distance_categories)))

    # Plotting
    grouped_data = data.groupby(['Probability Group', 'Distance Group']).size().unstack(fill_value=0)
    grouped_data.plot(kind='bar', stacked=True, figsize=(12, 8), color=colors)
    
    plt.title('Counts of Distance Groups Across Probability Percentages')
    plt.xlabel('Probability Percentage Groups')
    plt.ylabel('Count of Species')
    plt.legend(title='Distance Groups', bbox_to_anchor=(1.05, 1), loc='upper left', labels=[str(x) for x in distance_categories])
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()



plot_probability_histogram_by_distance_groups(path_nwk, path_tsv, species_ref)
# plot_distance_histogram_by_probability_groups(path_nwk, path_tsv, species_ref)



# plot_grouped_histogram_by_distance(path_nwk, path_tsv, species_ref)
# plot_distances(path_nwk, path_tsv, species_ref)
# plot_probability_curve_with_annotations(path_nwk, path_tsv, species_ref)
