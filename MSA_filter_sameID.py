from Bio import SeqIO

def remove_duplicates(msa_file):
    sequences = SeqIO.parse(msa_file, "fasta")
    unique_sequences = {}
    for seq in sequences:
        if seq.id not in unique_sequences:
            unique_sequences[seq.id] = seq
    return unique_sequences

def filter_and_match_msa(msa_file1, msa_file2, output_file1, output_file2, excluded_output_file1, excluded_output_file2):
    msa1_unique = remove_duplicates(msa_file1)
    msa2_unique = remove_duplicates(msa_file2)

    common_sequences1 = {}
    common_sequences2 = {}
    uniq_sequences1 = {k: v for k, v in msa1_unique.items() if k not in msa2_unique}
    uniq_sequences2 = {k: v for k, v in msa2_unique.items() if k not in msa1_unique}
    
    # Trouver les séquences communes
    for seq_id in msa1_unique:
        if seq_id in msa2_unique and not seq_id.startswith("EN"):
            common_sequences1[seq_id] = msa1_unique[seq_id]
            common_sequences2[seq_id] = msa2_unique[seq_id]

    # Écrire les séquences communes
    with open(output_file1, "w") as output_handle1, open(output_file2, "w") as output_handle2:
        SeqIO.write(common_sequences1.values(), output_handle1, "fasta")
        SeqIO.write(common_sequences2.values(), output_handle2, "fasta")

    # Écrire les séquences uniques (exclues)
    with open(excluded_output_file1, "w") as excluded_handle1, open(excluded_output_file2, "w") as excluded_handle2:
        SeqIO.write(uniq_sequences1.values(), excluded_handle1, "fasta")
        SeqIO.write(uniq_sequences2.values(), excluded_handle2, "fasta")
    
    print(f"Fichiers créés : '{output_file1}', '{output_file2}', '{excluded_output_file1}', et '{excluded_output_file2}' avec les séquences nécessaires.")







# Exemple d'utilisation
msa_file1 = "/Users/louiscarrel/Documents/Alignement_Project/largescale_kinase/DATA/ENSG00000010810/1st_transcrit_5_0/New_alignement/msa_s_exon_5_0.fasta"
msa_file2 = "/Users/louiscarrel/Documents/Alignement_Project/largescale_kinase/DATA/ENSG00000010810/1st_transcrit_5_0/New_alignement/msa_s_exon_5_1.fasta"
output_file1 = "/Users/louiscarrel/Documents/Alignement_Project/largescale_kinase/DATA/ENSG00000010810/1st_transcrit_5_0/New_alignement/undecided_can.fasta"
output_file2 = "/Users/louiscarrel/Documents/Alignement_Project/largescale_kinase/DATA/ENSG00000010810/1st_transcrit_5_0/New_alignement/undecided_alt.fasta"
excluded_output_file1 = "/Users/louiscarrel/Documents/Alignement_Project/largescale_kinase/DATA/ENSG00000010810/1st_transcrit_5_0/New_alignement/excluded_can.fasta"
excluded_output_file2 = "/Users/louiscarrel/Documents/Alignement_Project/largescale_kinase/DATA/ENSG00000010810/1st_transcrit_5_0/New_alignement/excluded_alt.fasta"
filter_and_match_msa(msa_file1, msa_file2, output_file1,output_file2,excluded_output_file1,excluded_output_file2)
