import sys
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

def merge_msas(msa_paths, output_file_all, output_file_ens):
    # Créer deux listes pour les séquences finales
    all_sequences = []
    ens_sequences = []

    # Parcourir tous les fichiers MSA donnés en entrée
    for path in msa_paths:
        # Lire le fichier MSA
        msa = AlignIO.read(path, "fasta")
        # Filtrer les séquences
        for record in msa:
            if record.id.startswith("ENS"):
                ens_sequences.append(record)
            else:
                all_sequences.append(record)

    # Créer les objets MSA à partir des listes de séquences
    all_msa = MultipleSeqAlignment(all_sequences)
    ens_msa = MultipleSeqAlignment(ens_sequences)

    # Écrire les nouveaux fichiers MSA
    AlignIO.write(all_msa, output_file_all, "fasta")
    AlignIO.write(ens_msa, output_file_ens, "fasta")

def main():
    # Lire les chemins des fichiers MSA et les fichiers de sortie de la ligne de commande
    if len(sys.argv) < 6:
        print("Usage: python script.py <msa1> <msa2> <msa3> <output_all.fasta> <output_ens.fasta>")
        sys.exit(1)

    msa_paths = sys.argv[1:4]
    output_file_all = sys.argv[4]
    output_file_ens = sys.argv[5]

    merge_msas(msa_paths, output_file_all, output_file_ens)

if __name__ == "__main__":
    main()
