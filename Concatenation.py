from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
import os
from Bio import SeqIO


base_path = "DATA/100species/ENSG00000010810/thoraxe/msa/"
file_names = [os.path.join(base_path, f"msa_s_exon_{i}.fasta") for i in ("38_18","38_19","38_20","38_21","38_22","38_23")]


# Charger les alignements
alignments = [AlignIO.read(filename, "fasta") for filename in file_names]

# Trouver les IDs communs à tous les fichiers
common_ids = set(record.id for record in alignments[0])
for align in alignments[1:]:
    common_ids.intersection_update(set(record.id for record in align))

# Filtrer chaque alignement pour ne garder que les séquences avec des IDs communs
filtered_alignments = []
for align in alignments:
    filtered = [record for record in align if record.id in common_ids]
    filtered_alignments.append(filtered)

# Créer un nouvel alignement avec les séquences concaténées
concatenated_records = []
for records in zip(*filtered_alignments):
    new_seq = "".join(str(record.seq) for record in records)  # Convertir Seq en str avant de concaténer
    concatenated_records.append(SeqRecord(Seq(new_seq), id=records[0].id, description=""))

# Créer un objet MultipleSeqAlignment à partir des SeqRecords
concatenated_alignment = MultipleSeqAlignment(concatenated_records)

# Enregistrer le nouvel alignement concaténé
AlignIO.write(concatenated_alignment, os.path.join("DATA/100species/ENSG00000010810/thoraxe/msa", "combined_filtered_alt.fasta"), "fasta")

# def remove_first_letter(input_file, output_file):
#     with open(input_file, "r") as infile, open(output_file, "w") as outfile:
#         for record in SeqIO.parse(infile, "fasta"):
#             modified_sequence = record.seq[1:]  # Supprime la première lettre
#             record.seq = modified_sequence
#             SeqIO.write(record, outfile, "fasta")


# # Liste des identifiants à conserver
# ids_to_keep = {"ENSBTAG00000007876", "ENSG00000107643", "ENSGGOG00000011771", "ENSMMUG00000004060", 
#                "ENSMODG00000002193", "ENSMUSG00000021936", "ENSOANG00000012095", "ENSRNOG00000020155",
#                "ENSSSCG00000010380", "ENSXETG00000021691"}

# remove_first_letter("ENSG00000107643/thoraxe_2/combined_filtered_alt.fasta","ENSG00000107643/thoraxe_2/combined_filtered_minus_1_alt.fasta")
# # Chemin des fichiers de sortie à lire
# file_paths = ["ENSG00000107643/thoraxe_2/combined_filtered_can.fasta", "ENSG00000107643/thoraxe_2/combined_filtered_minus_1_alt.fasta"]

# non_kept_sequences = []


# for file_path in file_paths:
#     # Lire l'alignement original
#     alignment = AlignIO.read(file_path, "fasta")
    
#     # Filtrer pour ne garder que les séquences avec les IDs spécifiés
#     filtered_alignment = [record for record in alignment if record.id in ids_to_keep]
    
#     # Créer un nouvel objet MultipleSeqAlignment avec les séquences filtrées
#     new_alignment = MultipleSeqAlignment(filtered_alignment)
    
#     # Construire le nom du fichier de sortie basé sur le fichier d'entrée
#     output_file = f"ENSG00000107643/thoraxe_2/filtered_{os.path.basename(file_path)}"
    
#     # Enregistrer le nouvel alignement dans un fichier
#     AlignIO.write(new_alignment, output_file, "fasta")



# for file_path in file_paths:
#     # Lire l'alignement original
#     alignment = AlignIO.read(file_path, "fasta")
    
#     # Collecter les séquences non conservées
#     non_kept = [record for record in alignment if record.id not in ids_to_keep]
#     non_kept_sequences.extend(non_kept)

# # Définir le nom du fichier de sortie pour les séquences non conservées
# output_file = "ENSG00000107643/thoraxe_2/non_kept_sequences.fasta"

# # Enregistrer les séquences non conservées dans un fichier FASTA
# with open(output_file, "w") as f:
#     for record in non_kept_sequences:
#         f.write(f">{record.id}\n{str(record.seq)}\n")





