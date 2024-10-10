"""
Script : Fasta_for_a3M.py
Description:
This script generates FASTA files from transcript data for a given gene. It is used before running ColabFold to prepare multiple sequence 
alignments in A3M format. The script reads two CSV files: path_table.csv and s_exon_table.csv, which contain information about exon paths and exon sequences, 
respectively. It then generates FASTA files, where each file corresponds to a transcript and contains the concatenated sequence of the exons.


"""
import os
import sys
import csv

def generate_fasta_for_transcripts(GENE, output_dir, transcript_ids=None):
    """
    Génère des fichiers FASTA pour les transcrits d'un gène donné. 
    Si des transcript_ids sont spécifiés, ne génère que ceux-là; sinon, génère tous les transcrits.
    Les fichiers FASTA seront créés dans un répertoire organisé par gène et transcrit.
    """
    path_table = []
    s_exon = {}

    print(f"Reading path_table.csv from {os.path.join(GENE, 'thoraxe/path_table.csv')}")

    # Lire path_table.csv
    try:
        with open(os.path.join(GENE, "thoraxe/path_table.csv"), 'r') as path_file:
            reader = csv.DictReader(path_file)
            for row in reader:
                path_table.append(row)
    except FileNotFoundError:
        print(f"Error: Could not find path_table.csv in {os.path.join(GENE, 'thoraxe')}")
        return

    print(f"Reading s_exon_table.csv from {os.path.join(GENE, 'thoraxe/s_exon_table.csv')}")
    # Lire s_exon_table.csv et construire un dictionnaire qui associe chaque exon spécifique à un gène et transcrit
    try:
        with open(os.path.join(GENE, "thoraxe/s_exon_table.csv"), 'r') as s_exon_file:
            reader = csv.DictReader(s_exon_file)
            for row in reader:
                key = (row['GeneID'], row['TranscriptIDCluster'], row['S_exonID'])  # Clé basée sur GeneID, TranscriptIDCluster et S_exonID
                s_exon[key] = row['S_exon_Sequence']  # Stocker la séquence de l'exon spécifique à ce gène et transcrit
    except FileNotFoundError:
        print(f"Error: Could not find s_exon_table.csv in {os.path.join(GENE, 'thoraxe')}")
        return

    # Itérer sur chaque ligne de path_table
    for row in path_table:
        gene_id = row['GeneID']
        transcript_id = row['TranscriptIDCluster'].replace('/', '_')  # Remplacer '/' par '_'

        # Si des transcript_ids sont spécifiés, ne traiter que ceux-là
        if transcript_ids and transcript_id not in transcript_ids:
            continue  # Passer ce transcript s'il n'est pas dans la liste

        # Créer des répertoires pour le gène et la transcription
        gene_dir = os.path.join(output_dir, gene_id)
        transcript_dir = os.path.join(gene_dir, transcript_id)

        print(f"Creating directory: {transcript_dir}")
        os.makedirs(transcript_dir, exist_ok=True)  # Crée le répertoire s'il n'existe pas

        fasta_filename = os.path.join(transcript_dir, 'transcript.fasta')

        # Créer un nouveau fichier FASTA dans le répertoire du transcrit
        print(f"Writing FASTA file for {transcript_id} to {fasta_filename}")
        with open(fasta_filename, 'w') as fasta_file:
            # Écrire l'en-tête FASTA pour ce transcrit
            fasta_file.write(f">{gene_id}-{transcript_id}\n")

            # Obtenir la liste des ID d'exon à partir de la colonne Path
            exon_list = row['Path'].split('/')
            exon_list = [exon for exon in exon_list if exon not in ['start', 'stop'] and not exon.startswith('0_')]

            # Concaténer les séquences des exons spécifiques au gène et transcrit
            final_sequence = ''
            for exon_id in exon_list:
                key = (gene_id, row['TranscriptIDCluster'], exon_id)  # Vérifier pour ce gène, ce transcrit et cet exon spécifique
                exon_sequence = s_exon.get(key)

                if exon_sequence:
                    final_sequence += exon_sequence
                else:
                    print(f"Warning: Séquence pour l'exon {exon_id} du gène {gene_id} et transcrit {transcript_id} introuvable.")
                    final_sequence += 'NNNN'  # Ajouter un placeholder si l'exon est manquant

            # Écrire la séquence concaténée dans le fichier FASTA
            fasta_file.write(final_sequence + '\n')

        info_filename = os.path.join(transcript_dir, 'info.txt')
        transcript_id_with_dashes = transcript_id.replace('_', '/')
        print(f"Writing info.txt for {transcript_id_with_dashes} to {info_filename}")

        with open(info_filename, 'w') as info_file:
            # Écrire l'ID du transcrit avec les tirets dans le fichier info.txt
            info_file.write(f"GeneID: {gene_id}\n")
            info_file.write(f"TranscriptID: {transcript_id_with_dashes}\n")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 Fasta_for_a3m.py <GENE_path> <output_dir> [TranscriptIDCluster1 TranscriptIDCluster2 ...]")
        sys.exit(1)

    GENE_path = sys.argv[1]  # Répertoire contenant les fichiers CSV
    output_dir = sys.argv[2]  # Répertoire racine de sortie

    # Récupérer les TranscriptIDCluster passés en arguments, s'ils existent
    if len(sys.argv) > 3:
        transcript_ids = [tid.replace('/', '_') for tid in sys.argv[3:]]
    else:
        transcript_ids = None  # Si aucun transcript n'est spécifié, générer tous les transcripts

    # Générer le fichier FASTA pour les transcrits
    print(f"Generating FASTA files for gene directory: {GENE_path}")
    generate_fasta_for_transcripts(GENE_path, output_dir, transcript_ids)
    print("FASTA file generation complete.")