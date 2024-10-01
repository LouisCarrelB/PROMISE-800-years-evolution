""" 
Permet de générer une liste de transcrit enregistré sous format texte afin  tout les chemis intéressant de l'ESG
"""
import pandas as pd
import sys

def compare_exons(canonical_path, alternative_path):
    # Diviser les chemins en segments d'exons
    canonical_exons = canonical_path.split('/')
    alternative_exons = alternative_path.split('/')
    
    # Initialiser des listes pour les exons différents
    different_canonical = []
    different_alternative = []
    
    # Comparer les exons, en tenant compte des longueurs différentes
    max_length = max(len(canonical_exons), len(alternative_exons))
    
    for i in range(max_length):
        # Gérer les cas où une liste est plus courte que l'autre
        canonical_exon = canonical_exons[i] if i < len(canonical_exons) else None
        alternative_exon = alternative_exons[i] if i < len(alternative_exons) else None
        
        if canonical_exon != alternative_exon:
            if canonical_exon:
                different_canonical.append(canonical_exon)
            if alternative_exon:
                different_alternative.append(alternative_exon)
    
    return different_canonical, different_alternative

def compare_genes(canonical_genes, alternative_genes):
    # Diviser les gènes en ID séparés par '/'
    canonical_genes_set = set(canonical_genes.split('/')) if canonical_genes else set()
    alternative_genes_set = set(alternative_genes.split('/')) if alternative_genes else set()
    
    # Comparer et trouver les gènes qui ne sont que dans CanonicalPath ou AlternativePath
    exclusive_to_canonical = canonical_genes_set - alternative_genes_set
    exclusive_to_alternative = alternative_genes_set - canonical_genes_set
    
    return exclusive_to_canonical, exclusive_to_alternative

def select_genes(exclusive_canonical_genes, exclusive_alternative_genes):
    # Convertir les ensembles en listes pour pouvoir accéder aux éléments
    exclusive_canonical_genes = list(exclusive_canonical_genes)
    exclusive_alternative_genes = list(exclusive_alternative_genes)

    # Initialiser les gènes décisifs pour CanonicalPath et AlternativePath
    decisive_canonical_gene = None
    decisive_alternative_gene = None

    # Si une des listes est vide, on prend ENSG pour la liste vide, sinon on prend le premier gène différent
    if not exclusive_canonical_genes:
        decisive_alternative_gene = next((gene for gene in exclusive_alternative_genes if gene.startswith("ENSG")), exclusive_alternative_genes[0])
        decisive_canonical_gene = "ENSG00000010810"
    elif not exclusive_alternative_genes:
        decisive_canonical_gene = next((gene for gene in exclusive_canonical_genes if gene.startswith("ENSG")), exclusive_canonical_genes[0])
        decisive_alternative_gene = "ENSG00000010810"
    
    # Si les deux listes sont remplies
    else:
        # On cherche d'abord un ENSG dans chaque liste
        decisive_canonical_gene = next((gene for gene in exclusive_canonical_genes if gene.startswith("ENSG")), exclusive_canonical_genes[0])
        decisive_alternative_gene = next((gene for gene in exclusive_alternative_genes if gene.startswith("ENSG")), exclusive_alternative_genes[0])
        
        # Si les deux gènes sont identiques, on choisit un gène différent pour l'Alternative
        if decisive_canonical_gene == decisive_alternative_gene:
            decisive_alternative_gene = exclusive_alternative_genes[0] if decisive_canonical_gene != exclusive_alternative_genes[0] else exclusive_alternative_genes[1]

    return decisive_canonical_gene, decisive_alternative_gene

def get_first_alternative_and_compare(path_csv, output_txt):
    # Charger le fichier CSV dans un DataFrame
    df = pd.read_csv(path_csv)
    
    # Chercher la première ligne où la colonne 'ASE' contient "alternative"
    first_alternative_row = df[df['ASE'].str.contains("alternative", na=False)].head(1)
    
    if first_alternative_row.empty:
        with open(output_txt, 'w') as file:
            file.write("Aucune ligne ne contient 'alternative' dans la colonne ASE.\n")
        return
    
    # Extraire les valeurs des colonnes CanonicalPath et AlternativePath
    canonical_path = first_alternative_row['CanonicalPath'].values[0]
    alternative_path = first_alternative_row['AlternativePath'].values[0]
    
    # Extraire les gènes de CanonicalPathGenes et AlternativePathGenes
    canonical_genes = first_alternative_row['CanonicalPathGenes'].values[0]
    alternative_genes = first_alternative_row['AlternativePathGenes'].values[0] if 'AlternativePathGenes' in first_alternative_row else ""

    # Comparer les exons
    different_canonical, different_alternative = compare_exons(canonical_path, alternative_path)
    
    # Comparer les gènes exclusifs
    exclusive_canonical_genes, exclusive_alternative_genes = compare_genes(canonical_genes, alternative_genes)
    
    # Sélectionner les gènes décisifs pour CAN et ALT
    decisive_canonical_gene, decisive_alternative_gene = select_genes(exclusive_canonical_genes, exclusive_alternative_genes)

    # Écrire les chemins, différences d'exons, gènes exclusifs et gènes sélectionnés dans un fichier texte
    with open(output_txt, 'w') as file:
        file.write(f"CanonicalPath: {canonical_path}\n")
        file.write(f"AlternativePath: {alternative_path}\n\n")
        
        if different_canonical or different_alternative:
            file.write(f"CAN : {'/'.join(different_canonical)}\n")
            file.write(f"ALT : {'/'.join(different_alternative)}\n")
        else:
            file.write("Aucune différence entre les exons.\n")
        
        # Enregistrer les gènes exclusifs avant de sélectionner les gènes décisifs
        file.write("\nGènes exclusifs CanonicalPath: " + '/'.join(exclusive_canonical_genes) + '\n')
        file.write("Gènes exclusifs AlternativePath: " + '/'.join(exclusive_alternative_genes) + '\n')
        
        # Écrire les gènes sélectionnés
        if decisive_canonical_gene and decisive_alternative_gene:
            file.write(f"\nCanonicalPathGenes: {decisive_canonical_gene}\n")
            file.write(f"AlternativePathGenes: {decisive_alternative_gene}\n")
        else:
            file.write("\nAucun gène sélectionné.\n")
    
    print(f"Comparaison terminée. Résultats enregistrés dans {output_txt}.")

def find_transcripts_by_paths(path_csv, ases_txt):
    # Charger le fichier path_table.csv dans un DataFrame
    df_path = pd.read_csv(path_csv)
    
    # Charger les informations depuis ases.txt
    with open(ases_txt, 'r') as file:
        lines = file.readlines()
    
    # Extraire les informations CanonicalPath et AlternativePath
    canonical_path = None
    alternative_path = None
    canonical_gene = None
    alternative_gene = None
    
    for line in lines:
        if line.startswith("CanonicalPath:"):
            canonical_path = line.split(":")[1].strip()
        elif line.startswith("AlternativePath:"):
            alternative_path = line.split(":")[1].strip()
        elif line.startswith("CanonicalPathGenes:"):
            canonical_gene = line.split(":")[1].strip()
        elif line.startswith("AlternativePathGenes:"):
            alternative_gene = line.split(":")[1].strip()

    if not canonical_path or not alternative_path or not canonical_gene or not alternative_gene:
        print("Erreur : Impossible d'extraire les chemins ou les gènes du fichier ases.txt.")
        return
    
    # Rechercher les transcrits correspondant au CanonicalPath et AlternativePath
    canonical_transcript = df_path[(df_path['GeneID'] == canonical_gene) & (df_path['Path'].str.contains(canonical_path))]
    alternative_transcript = df_path[(df_path['GeneID'] == alternative_gene) & (df_path['Path'].str.contains(alternative_path))]

    # Extraire les ID des transcrits
    transcriptID_CAN = canonical_transcript['TranscriptIDCluster'].values[0] if not canonical_transcript.empty else "Aucun"
    transcriptID_ALT = alternative_transcript['TranscriptIDCluster'].values[0] if not alternative_transcript.empty else "Aucun"

    # Ouvrir le fichier ases.txt en mode écriture ("w" écrase tout)
    with open(ases_txt, 'w') as file:
        # Réécrire le contenu original
        file.writelines(lines)
        
        # Ajouter les transcritID_CAN et transcritID_ALT à la fin du fichier
        file.write("\ntranscritID_CAN : " + transcriptID_CAN + "\n")
        file.write("transcritID_ALT : " + transcriptID_ALT + "\n")

    print(f"Les résultats ont été ajoutés et écrits dans {ases_txt}.")

if __name__ == '__main__':

    path_csv = sys.argv[1] + "/thoraxe/ases_table.csv"
    path_path = sys.argv[1] + "/thoraxe/path_table.csv"
    output_txt = sys.argv[2] + "/ases.txt"
    get_first_alternative_and_compare(path_csv, output_txt)
    find_transcripts_by_paths(path_path, output_txt)