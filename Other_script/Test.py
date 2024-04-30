import pandas as pd 

def find_and_update_gene_ids(dataframe, ases_table_path):
    # Charger la dataframe ases_table et lire uniquement la première ligne
    ases_df = pd.read_csv(ases_table_path, nrows=1)
    
    # Lire les chemins depuis les colonnes
    alt_path = ases_df['AlternativePath'].iloc[0]
    can_path = ases_df['CanonicalPath'].iloc[0]
    
    # Fonction pour extraire les GeneID pour un chemin donné
    def find_gene_ids_with_path(s_exon_path):
        path_segments = s_exon_path.split('/')
        valid_gene_ids = set(dataframe['GeneID'])
        
        for segment in path_segments:
            current_segment_gene_ids = set(dataframe[dataframe['S_exonID'] == segment]['GeneID'])
            valid_gene_ids.intersection_update(current_segment_gene_ids)
        
        # Filtrer les GeneID pour exclure ceux commençant par 'ENS'
        filtered_ids = [gene_id for gene_id in valid_gene_ids if not gene_id.startswith('ENS')]
        return '/'.join(filtered_ids)  # Convertir la liste en chaîne formatée
    
    # Récupérer les GeneID pour chaque chemin
    alt_gene_ids = find_gene_ids_with_path(alt_path)
    can_gene_ids = find_gene_ids_with_path(can_path)
    
    # Calculer les GeneID communs aux deux listes
    common_gene_ids = '/'.join(set(alt_gene_ids.split('/')) & set(can_gene_ids.split('/')))
    
    # Ajouter les résultats aux colonnes spécifiques de ases_df
    ases_df['CanonicalPathGenes'] = ases_df.get('CanonicalPathGenes', '') + can_gene_ids
    ases_df['AlternativePathGenes'] = ases_df.get('AlternativePathGenes', '') + alt_gene_ids
    ases_df['CommonGenes'] = ases_df.get('CommonGenes', '') + common_gene_ids
    
    # Sauvegarder les modifications dans un nouveau fichier
    output_path = ases_table_path.replace('ases_table.csv', 'ases_table_a3m.csv')
    ases_df.to_csv(output_path, index=False)

    return alt_gene_ids, can_gene_ids, common_gene_ids

# Charger les données de s_exon_table
dataframe = pd.read_csv('../DATA/ENSG00000107643/thoraxe/s_exon_table_a3m.csv')

# Appeler la fonction et imprimer les résultats
alt_ids, can_ids, common_ids = find_and_update_gene_ids(dataframe, '../DATA/ENSG00000107643/thoraxe/ases_table.csv')
print("Alt Gene IDs:", alt_ids)
print("Can Gene IDs:", can_ids)
print("Common Gene IDs:", common_ids)
