import os
from Bio import SeqIO
from collections import defaultdict
import sys

def combine_msas(input_dir):
    # Créer le répertoire "all_msas/" dans le même répertoire que l'input_dir
    output_dir = os.path.join(input_dir, "all_msas")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Dictionnaire pour stocker les séquences par nom de fichier
    msa_sequences = defaultdict(lambda: defaultdict(set))
    
    # Parcourir les fichiers d'entrée pour traiter les MSAs
    for subdir, _, files in os.walk(input_dir):
        for file in files:
            if file.endswith(".fasta") or file.endswith(".msa"):  
                filepath = os.path.join(subdir, file)
                
                # Lire et stocker les séquences pour chaque fichier
                with open(filepath, "r") as f:
                    for record in SeqIO.parse(f, "fasta"):
                        msa_sequences[file][str(record.seq)].add(record.description)
    
    # Écriture des MSAs combinés et dédupliqués dans "all_msas/"
    for msa_name, sequences in msa_sequences.items():
        output_path = os.path.join(output_dir, msa_name)
        with open(output_path, "w") as output_file:
            for seq, descriptions in sequences.items():
                description = ";".join(descriptions)  # Concaténer les descriptions si nécessaire
                output_file.write(f">{description}\n{seq}\n")

    print(f"MSAs combinés et dédupliqués écrits dans : {output_dir}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 combine_msas.py <gene_name>")
        sys.exit(1)

    gene_name = sys.argv[1]
    gene_dir = os.path.join('/shared/home/carrell', gene_name)
    input_dir = os.path.join(gene_dir, "New_alignement/")
    
    combine_msas(input_dir)