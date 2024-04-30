import pandas as pd

def build_species_dict(exon_file):
    # Lire les données du fichier s_exon_table.csv
    exon_df = pd.read_csv(exon_file)
    
    # Créer un dictionnaire à partir de la colonne GeneID et Species, en formatant les noms d'espèces
    species_dict = pd.Series(
        exon_df['Species'].apply(lambda x: x.lower().capitalize()).values, # Met la première lettre de chaque mot en majuscule
        index=exon_df['GeneID']
    ).to_dict()
    return species_dict



def read_gene_data(file_path, row_number, species_dict):
    
    df = pd.read_csv(file_path)
    
    selected_row = df.iloc[row_number]
    
    canonical_genes = set(selected_row['CanonicalPathGenes'].split('/'))
    alternative_genes = set(selected_row['AlternativePathGenes'].split('/'))
    
    gene_dict = {gene: {'species': species_dict.get(gene, 'Unknown species'), 'index': 'Not_in_path'} for gene in species_dict.keys()}
    
    for gene in canonical_genes:
        if gene in gene_dict:
            gene_dict[gene]['index'] = 'can'
        else:
            gene_dict[gene] = {'species': species_dict.get(gene, 'Unknown species'), 'index': 'can'}
    
    for gene in alternative_genes:
        if gene in gene_dict:
            if gene_dict[gene]['index'] == 'can':
                gene_dict[gene]['index'] = 'both'
            else:
                gene_dict[gene]['index'] = 'alt'
        else:
            gene_dict[gene] = {'species': species_dict.get(gene, 'Unknown species'), 'index': 'alt'}
    
    result_df = pd.DataFrame.from_dict(gene_dict, orient='index', columns=['species', 'index'])
    
    # Assurer une ligne par espèce
    result_df = result_df.drop_duplicates(subset='species').reset_index(drop=True)
    
    return result_df

# Chemins des fichiers - Assurez-vous que les chemins sont corrects
file_path = 'ENSG00000107643/thoraxe_2/ases_table.csv'
exon_file = 'ENSG00000107643/thoraxe_2/s_exon_table.csv'

# Construire le dictionnaire d'espèces
species_dict = build_species_dict(exon_file)

# Utiliser la fonction
row_number = 0  
output_table = read_gene_data(file_path, row_number, species_dict)
print(output_table)

output_table.to_csv('ENSG00000107643/thoraxe_2/species_path_table.tsv', sep='\t', index=False, columns=['species', 'index'])