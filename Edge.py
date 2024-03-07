from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pds
from io import StringIO
from Bio import pairwise2
import sys



def filter_sequences(a3m, output_file):
    with open(a3m, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            if line.startswith('>'):  
                f_out.write(line)
            else:
                filtered_sequence = ''.join(['-' if c != '-' and not c.islower() else c for c in line.strip()])
                f_out.write(filtered_sequence + '\n')


def filter_sequences_with_positions(input_file, output_file):
    sequences_with_positions = []

    with open(input_file, 'r') as f_in:
        current_sequence = ""
        current_start = 0

        for line in f_in:
            if line.startswith('>'):  
                if current_sequence:
                    sequences_with_positions.append((current_sequence, current_start, len(current_sequence)))
                current_sequence = ""  # Réinitialiser la séquence actuelle
                continue

            filtered_sequence = ''.join(['-' if c != '-' and not c.islower() else c for c in line.strip()])
            for i, char in enumerate(filtered_sequence):
                if char.islower():  # Si le caractère est une lettre minuscule
                    if not current_sequence:
                        current_start = i
                    current_sequence += char
                elif current_sequence:  # Si le caractère n'est pas une lettre minuscule mais la séquence actuelle n'est pas vide
                    sequences_with_positions.append((current_sequence, current_start, i))
                    current_sequence = ""  # Réinitialiser la séquence actuelle

        # Ajouter la dernière séquence
        if current_sequence:
            sequences_with_positions.append((current_sequence, current_start, len(current_sequence)))

    # Écrire les séquences filtrées avec les positions dans le fichier de sortie
    with open(output_file, 'w') as f_out:
        for sequence, start, end in sequences_with_positions:
            f_out.write(f'{sequence} [{start},{end}]\n')


def filter_sequences_with_positions_and_ids(input_file, output_file):
    sequences_with_positions_and_ids = []
    sequence_id = 0

    with open(input_file, 'r') as f_in:
        current_sequence = ""
        current_start = 0
        line_number = 0

        for line in f_in:
            line_number += 1
            if line.startswith('>'): 
                if current_sequence:
                    sequences_with_positions_and_ids.append((current_sequence, current_start, len(current_sequence), start_line, sequence_id))
                    sequence_id += 1
                current_sequence = ""  # Réinitialiser la séquence actuelle
                continue

            filtered_sequence = ''.join(['-' if c != '-' and not c.islower() else c for c in line.strip()])
            for i, char in enumerate(filtered_sequence):
                if char.islower():  # Si le caractère est une lettre minuscule
                    if not current_sequence:
                        current_start = i
                        start_line = line_number
                    current_sequence += char
                elif current_sequence:  # Si le caractère n'est pas une lettre minuscule mais la séquence actuelle n'est pas vide
                    sequences_with_positions_and_ids.append((current_sequence, current_start, i, start_line, sequence_id))
                    current_sequence = ""  # Réinitialiser la séquence actuelle

        # Ajouter la dernière séquence
        if current_sequence:
            sequences_with_positions_and_ids.append((current_sequence, current_start, len(current_sequence), start_line, sequence_id))

    # Écrire les séquences filtrées avec les positions et les IDs dans le fichier de sortie
    with open(output_file, 'w') as f_out:
        for sequence, start, end, line, seq_id in sequences_with_positions_and_ids:
            f_out.write(f'Sequence ID: {seq_id} - {sequence} [{start},{end}] (found in line {line})\n')




a3m = "../DATA/ENSG00000107643/15333.a3m"
output_file_path = '../DATA/ENSG00000107643/inter/reduced.a3m'
output2_file_path = '../DATA/ENSG00000107643/inter/info_a3m_reduced.a3m'
filter_sequences(a3m, output_file_path)
filter_sequences_with_positions(a3m, output2_file_path)
filter_sequences_with_positions_and_ids(a3m, output2_file_path)