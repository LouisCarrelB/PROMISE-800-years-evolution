#!/usr/bin/env python3

import os
import sys
import subprocess
import shutil 
from Bio import SeqIO


def realign_sequences_with_mafft(input_fasta, output_fasta):
    """
    Réaligne les séquences en utilisant MAFFT.
    
    Args:
        input_fasta (str): Chemin vers le fichier FASTA à réaligner.
        output_fasta (str): Chemin vers le fichier de sortie où écrire l'alignement réaligné.
    """
    print(f"Réalignement des séquences dans {input_fasta} avec MAFFT.")
    mafft_cmd = ["mafft", "--auto", input_fasta]
    with open(output_fasta, "w") as output_handle:
        result = subprocess.run(mafft_cmd, stdout=output_handle, stderr=subprocess.PIPE, text=True)
    
    if result.returncode != 0:
        print(f"Erreur lors du réalignement des séquences : {result.stderr}")
    else:
        print(f"Alignement terminé avec succès et enregistré dans {output_fasta}")


def combine_fasta_files(can_fasta, alt_fasta, combined_fasta):
    """
    Combine les fichiers FASTA de CAN et ALT en un seul fichier FASTA pour alignement.
    
    Args:
        can_fasta (str): Chemin du fichier CAN_ori.fasta.
        alt_fasta (str): Chemin du fichier ALT_ori.fasta.
        combined_fasta (str): Chemin du fichier FASTA combiné.
    """
    with open(combined_fasta, "w") as outfile:
        # Ajouter les séquences de CAN
        for record in SeqIO.parse(can_fasta, "fasta"):
            SeqIO.write(record, outfile, "fasta")
        
        # Ajouter les séquences de ALT
        for record in SeqIO.parse(alt_fasta, "fasta"):
            SeqIO.write(record, outfile, "fasta")

# Le reste du script reste inchangé
def count_sequences(fasta_file):
    """
    Compte le nombre de séquences dans un fichier FASTA.
    
    Args:
        fasta_file (str): Chemin vers le fichier FASTA.
    
    Returns:
        int: Nombre total de séquences dans le fichier.
    """
    return sum(1 for _ in SeqIO.parse(fasta_file, "fasta"))

def split_combined_fasta_by_count(combined_fasta, can_fasta, alt_fasta, can_count):
    """
    Sépare les séquences alignées du fichier combiné en fichiers CAN et ALT en fonction du nombre de séquences.
    
    Args:
        combined_fasta (str): Chemin du fichier FASTA combiné aligné.
        can_fasta (str): Chemin du fichier CAN_ori_aligned.fasta.
        alt_fasta (str): Chemin du fichier ALT_ori_aligned.fasta.
        can_count (int): Nombre de séquences qui appartiennent à CAN.
    """
    with open(can_fasta, "w") as can_outfile, open(alt_fasta, "w") as alt_outfile:
        for i, record in enumerate(SeqIO.parse(combined_fasta, "fasta")):
            if i < can_count:  # Les premières `can_count` séquences vont dans CAN
                SeqIO.write(record, can_outfile, "fasta")
            else:  # Le reste va dans ALT
                SeqIO.write(record, alt_outfile, "fasta")

def create_seed_and_unified(can_aligned, alt_aligned, unified_fasta, seed_can_fasta, seed_alt_fasta):
    """
    Crée les fichiers seed_CAN.fasta, seed_ALT.fasta et unified.fasta.
    
    Args:
        can_aligned (str): Chemin vers le fichier aligné CAN_ori_aligned.fasta.
        alt_aligned (str): Chemin vers le fichier aligné ALT_ori_aligned.fasta.
        unified_fasta (str): Chemin vers le fichier unified.fasta à créer.
        seed_can_fasta (str): Chemin vers le fichier seed_CAN.fasta à créer.
        seed_alt_fasta (str): Chemin vers le fichier seed_ALT.fasta à créer.
    """
    seq_dict = {}  # Pour éviter les doublons dans unified
    seed_can_records = []
    seed_alt_records = []
    
    # Unifier et créer les fichiers seeds
    with open(unified_fasta, "w") as unified_handle:
        # Process CAN
        for record in SeqIO.parse(can_aligned, "fasta"):
            if record.id.startswith("ENS"):  # Ajouter aux seed_CAN
                seed_can_records.append(record)
            else:  # Ajouter à unified s'il n'est pas déjà présent
                seq_key = (record.id, str(record.seq))  # Utiliser ID et séquence comme clé pour éviter les doublons
                if seq_key not in seq_dict:
                    seq_dict[seq_key] = record
                    SeqIO.write(record, unified_handle, "fasta")

        # Process ALT
        for record in SeqIO.parse(alt_aligned, "fasta"):
            if record.id.startswith("ENS"):  # Ajouter aux seed_ALT
                seed_alt_records.append(record)
            else:  # Ajouter à unified s'il n'est pas déjà présent
                seq_key = (record.id, str(record.seq))
                if seq_key not in seq_dict:
                    seq_dict[seq_key] = record
                    SeqIO.write(record, unified_handle, "fasta")

    # Écriture des fichiers seed_CAN et seed_ALT
    with open(seed_can_fasta, "w") as seed_can_handle:
        SeqIO.write(seed_can_records, seed_can_handle, "fasta")
    with open(seed_alt_fasta, "w") as seed_alt_handle:
        SeqIO.write(seed_alt_records, seed_alt_handle, "fasta")

    print(f"Unified FASTA créé : {unified_fasta}")
    print(f"Seed CAN FASTA créé : {seed_can_fasta}")
    print(f"Seed ALT FASTA créé : {seed_alt_fasta}")

def unify_fasta_files(gene_dir):
    # Définir les chemins des fichiers
 
    base_path = os.path.join(gene_dir, "New_alignement/exclusif_paths/")
    can_fasta_path = os.path.join(base_path, "CAN_ori.fasta")
    alt_fasta_path = os.path.join(base_path, "ALT_ori.fasta")
    combined_fasta_path = os.path.join(base_path, "combined.fasta")
    combined_aligned_fasta_path = os.path.join(base_path, "combined_aligned.fasta")
    can_aligned_fasta_path = os.path.join(base_path, "CAN_ori_aligned.fasta")
    alt_aligned_fasta_path = os.path.join(base_path, "ALT_ori_aligned.fasta")
    unified_fasta_path = os.path.join(base_path, "unified.fasta")
    seed_can_fasta_path = os.path.join(base_path, "seed_CAN.fasta")
    seed_alt_fasta_path = os.path.join(base_path, "seed_ALT.fasta")

    # Compter les séquences pour CAN et ALT
    can_count = count_sequences(can_fasta_path)
    alt_count = count_sequences(alt_fasta_path)
    
    print(f"Nombre de séquences dans CAN : {can_count}")
    print(f"Nombre de séquences dans ALT : {alt_count}")

    # Étape 1 : Combiner les fichiers CAN et ALT en un seul fichier
    combine_fasta_files(can_fasta_path, alt_fasta_path, combined_fasta_path)
    print(f"Fichiers CAN et ALT combinés dans {combined_fasta_path}")

    # Étape 2 : Réaligner les séquences combinées avec MAFFT
    realign_sequences_with_mafft(combined_fasta_path, combined_aligned_fasta_path)

    # Étape 3 : Séparer le fichier combiné aligné en fichiers CAN et ALT alignés par compte
    split_combined_fasta_by_count(combined_aligned_fasta_path, can_aligned_fasta_path, alt_aligned_fasta_path, can_count)
    print(f"Fichiers alignés séparés : {can_aligned_fasta_path} et {alt_aligned_fasta_path}")

    # Étape 4 : Créer unified.fasta, seed_CAN.fasta et seed_ALT.fasta à partir des fichiers alignés
    create_seed_and_unified(can_aligned_fasta_path, alt_aligned_fasta_path, unified_fasta_path, seed_can_fasta_path, seed_alt_fasta_path)

    return unified_fasta_path,seed_can_fasta_path,seed_alt_fasta_path

def main(gene_dir,base_dir):
    # Unifier les fichiers FASTA (CAN et ALT)
    uni,can,alt = unify_fasta_files(gene_dir)
    print("process redistribution with weights")
    Redis_cmd = [
        'python3',
        os.path.join(base_dir, 'Alignement_2_msa.py'),
        uni,
        can,
        alt,
        os.path.join(gene_dir,'New_alignement','exclusif_paths')
    ]

    result = subprocess.run(Redis_cmd, capture_output=True, text=True)

    # Vérifier s'il y a eu une erreur lors de l'exécution
    if result.returncode != 0:
        print(result)
        print(result.stderr)
        print(result.stdout)  
    else:
        print(result.stdout) 

if __name__ == '__main__':
    if len(sys.argv) > 3:
        print("Usage: python3 workflow_redistribution.py <gene_name>")
        sys.exit(1)
    
    gene_name = sys.argv[1]
    gene_dir = os.path.join("/scratch/carrelbl/DATA/", gene_name)
    base_dir = sys.argv[2]
    main(gene_dir,base_dir)

  