import sys
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord

def trim_last_character(msa_path, output_path):
    # Charger le MSA original
    original_msa = AlignIO.read(msa_path, "fasta")
    
    # Créer un nouveau MSA sans la dernière lettre de chaque séquence
    trimmed_msa = MultipleSeqAlignment(
        SeqRecord(record.seq[:-1], id=record.id, description=record.description) 
        for record in original_msa
    )
    
    # Écrire le nouveau MSA dans un fichier
    AlignIO.write(trimmed_msa, output_path, "fasta")

def main():
    # Vérifier le nombre d'arguments
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_msa.fasta> <output_msa.fasta>")
        sys.exit(1)
    
    input_path = sys.argv[1]
    output_path = sys.argv[2]
    
    trim_last_character(input_path, output_path)

if __name__ == "__main__":
    main()
