"""
Script : Fasta_for_a3M.py
Description :
    Ce script génère des fichiers FASTA à partir des données de transcrits pour un gène donné.
    Il est utilisé avant l'exécution de ColabFold pour préparer les alignements multiples au
    format A3M. Le script lit deux fichiers CSV : 'path_table.csv' et 's_exon_table.csv',
    qui contiennent des informations sur les chemins d'exons et les séquences d'exons
    respectivement. Il génère ensuite des fichiers FASTA, où chaque fichier correspond à un
    transcrit et contient la séquence concaténée des exons.

    Fonctionnement :
    - Pour chaque gène et transcrit, un dossier est créé dans le répertoire 'results'.
    - Le script itère sur les exons d'un transcrit, concatène leurs séquences, et écrit le
      tout dans un fichier FASTA.
    - En cas de séquence manquante ou d'exon introuvable, le script ajoute une séquence
      placeholder ('NNNN') et signale le problème.

Utilisation :
    python Fasta_for_a3m.py <gene_name>

    Où <gene_name> correspond au nom du gène pour lequel les FASTA doivent être générés.

Exemple :
    python3 Fasta_for_a3m.py ENSG00000107643
"""


import os
import sys
import csv

def generate_fasta_for_transcripts(GENE, output_dir):
    # Charger les tables path_table et s_exon sans Pandas
    path_table = []
    s_exon = {}

    # Lire path_table.csv
    with open(os.path.join(GENE, "thoraxe/path_table.csv"), 'r') as path_file:
        reader = csv.DictReader(path_file)
        for row in reader:
            path_table.append(row)
    
    # Lire s_exon_table.csv
    with open(os.path.join(GENE, "thoraxe/s_exon_table.csv"), 'r') as s_exon_file:
        reader = csv.DictReader(s_exon_file)
        for row in reader:
            s_exon[row['S_exonID']] = row['S_exon_Sequence']

    # Itérer sur chaque ligne de path_table
    for row in path_table:
        gene_id = row['GeneID']
        transcript_id = row['TranscriptIDCluster'].replace('/', '_')  # Remplacer '/' par '_'
        
        # Créer des répertoires pour le gène et la transcription
        gene_dir = os.path.join(output_dir, gene_id)
        transcript_dir = os.path.join(gene_dir, transcript_id)
        os.makedirs(transcript_dir, exist_ok=True)  # Crée le répertoire s'il n'existe pas
        
        fasta_filename = os.path.join(transcript_dir, 'transcript.fasta')
        
        # Créer un nouveau fichier FASTA
        with open(fasta_filename, 'w') as fasta_file:
            # Écrire l'en-tête FASTA
            fasta_file.write(f">{gene_id}-{transcript_id}\n")
            
            # Obtenir la liste des ID d'exon à partir de la colonne Path
            exon_list = row['Path'].split('/')
            exon_list = [exon for exon in exon_list if exon not in ['start', 'stop'] and not exon.startswith('0_')]
            
            # Initialiser la séquence à écrire dans le fichier FASTA
            final_sequence = ''
            
            # Boucler sur chaque ID d'exon dans la liste
            for i, exon_id in enumerate(exon_list):
                exon_sequence = s_exon.get(exon_id)
                
                if exon_sequence:
                    final_sequence += exon_sequence
                else:
                    # Ajouter un placeholder si l'exon est manquant
                    print(f"Exon {exon_id} not found at position {i+1} in {gene_id}-{transcript_id}")
                    final_sequence += 'NNNN'
            
            # Écrire la séquence concaténée dans le fichier FASTA
            fasta_file.write(final_sequence + '\n')
        
        # Créer un fichier texte 'info' contenant l'ID du FASTA
        info_filename = os.path.join(transcript_dir, 'info.txt')
        with open(info_filename, 'w') as info_file:
            info_file.write(f"{gene_id}-{transcript_id}\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 Fasta_for_a3m.py <GENE_path> <output_dir>")
        sys.exit(1)
    
    GENE_path = sys.argv[1]
    output_dir = sys.argv[2]
    generate_fasta_for_transcripts(GENE_path, output_dir)