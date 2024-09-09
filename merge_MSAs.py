import os
from Bio import SeqIO
from collections import defaultdict

def combine_msas(input_dir, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Dictionnaire pour stocker les séquences par nom de fichier
    msa_sequences = defaultdict(lambda: defaultdict(set))
    
    
    for subdir, _, files in os.walk(input_dir):
        for file in files:
            if file.endswith(".fasta") or file.endswith(".msa"):  
                filepath = os.path.join(subdir, file)
                
                
                with open(filepath, "r") as f:
                    for record in SeqIO.parse(f, "fasta"):
                        msa_sequences[file][str(record.seq)].add(record.description)
    
    # Écriture des MSAs combinés et dédupliqués dans le répertoire de sortie
    for msa_name, sequences in msa_sequences.items():
        output_path = os.path.join(output_dir, msa_name)
        with open(output_path, "w") as output_file:
            for seq, descriptions in sequences.items():
                description = ";".join(descriptions)  # Concaténer les descriptions si nécessaire
                output_file.write(f">{description}\n{seq}\n")

    print(f"MSAs combinés et dédupliqués écrits dans : {output_dir}")

if __name__ == "__main__":
    input_directory = "MSAs/"
    output_directory = "New_alignementv2"
    combine_msas(input_directory, output_directory)