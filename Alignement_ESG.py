from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import glob
import re 
import os
import subprocess
from io import StringIO
from Bio import pairwise2
import sys
import shutil
import math 
import csv 
import tempfile
from tqdm import tqdm 


def write_msa_to_temp_file(msa_object, temp_file_path):
    """
    Writes the MSA object to a temporary file.

    Parameters:
    msa_object (list): A list of SeqRecord objects representing the MSA.
    temp_file_path (str): Path to the temporary file to write the MSA object.
    """
    SeqIO.write(msa_object, temp_file_path, "fasta")

def build_hmm(input_file, output_file):
    hmmbuild_command = f"/Users/louiscarrel/Downloads/hmmer-3.4/src/hmmbuild {output_file} {input_file}"
    subprocess.call(hmmbuild_command, shell=True,stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def score_sequence_with_hmm(hmm_file, sequence,exon_id):
    sequence_file = nouveau_repertoire +  f"temp_sequence_{exon_id}.fasta"
    # Remove dashes from the sequence
    cleaned_sequence = sequence.replace("-", "")

    # Save the new sequence in a FASTA file
    with open(sequence_file, 'w') as f:
        f.write(">new_sequence\n" + cleaned_sequence + "\n")

    score_output_path = os.path.join(nouveau_repertoire, f"score_output_{exon_id}.txt")
    # Execute hmmsearch
    hmmsearch_command = f"/Users/louiscarrel/Downloads/hmmer-3.4/src/hmmsearch --tblout {score_output_path} {hmm_file} {sequence_file}"
    subprocess.call(hmmsearch_command, shell=True,stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    os.remove(sequence_file)
    # Read and return the score from score_output.txt
    try:
        with open(nouveau_repertoire + f"score_output_{exon_id}.txt", 'r') as f:
            for line in f:
                if not line.startswith("#"):
                    return float(line.split()[5])  # needed to transfrom string in float 
    except FileNotFoundError:
        print("Score output file not found.")
        return None




def traiter_fichier_a3m(chemin_fichier):
    '''
    Traite un fichier au format A3M pour compter les occurrences des identifiants de séquence et ajuster les noms de transcrits en cas de doublons.
    
    Parameters:
    chemin_fichier (str): Chemin d'accès au fichier A3M à traiter.
    
    Returns:
    None: Cette fonction génère un fichier "good.a3m" avec des identifiants de séquence ajustés si nécessaire.
    '''
    occurrences = {}
    
    # Compter les occurrences des identifiants dans le fichier A3M
    with open(chemin_fichier, 'r') as file:
        for line in file:
            if line.startswith(">"):
                identifiant = line.split()[0][1:]
                if identifiant in occurrences:
                    occurrences[identifiant] += 1
                else:
                    occurrences[identifiant] = 1
    
    noms_transcrits = {}
    
    # Générer un nouveau fichier A3M avec des identifiants ajustés pour les doublons
    with open(chemin_fichier, 'r') as input_file, open(GENE + "good.a3m", 'w') as output_file:
        for line in input_file:
            if line.startswith(">"):
                identifiant = line.split()[0][1:]
                if occurrences[identifiant] > 1:
                    if identifiant not in noms_transcrits:
                        noms_transcrits[identifiant] = 1
                    else:
                        noms_transcrits[identifiant] += 1
                    nouveau_nom = f"{identifiant}.{noms_transcrits[identifiant]}"
                    output_file.write(">" + nouveau_nom + "\n")
                else:
                    output_file.write(line)
            else:
                output_file.write(line)

def gap_inter(A, B):
    '''
    Insère des tirets dans la séquence B aux positions correspondant aux tirets présents dans la séquence A.
    
    Parameters:
    A (str): La séquence de référence contenant les tirets.
    B (str): La séquence dans laquelle les tirets doivent être insérés.
    
    Returns:
    str: La séquence B ajustée avec des tirets insérés aux positions correspondantes à ceux de la séquence A.
    '''
    index_tiret = [i+1 for i in range(len(A)-2) if A[i].isalpha() and A[i+1] == '-' and A[i+2].isalpha()]
    for pos in reversed(index_tiret):
        if pos < len(B):
            B = B[:pos] + '-' + B[pos:]
    return B

def get_sexon_coord(gid, tid, seqFname):
    '''
    Extrait les coordonnées des exons à partir du fichier de séquence au format PIR en se basant sur l'identifiant du gène et du transcrit.
    
    Parameters:
    gid (str): Identifiant du gène.
    tid (str): Identifiant du transcrit.
    seqFname (str): Chemin d'accès au fichier de séquence PIR.
    
    Returns:
    list: Une liste de tuples représentant les coordonnées des exons (type d'exon, position de début, position de fin).
    '''
    fseq = open(seqFname)
    lines = fseq.readlines()
    fseq.close()

    i = 0
    found = False
    while (i < len(lines)) and (not found):
        found = lines[i].startswith('>P1;'+gid+' '+tid)
        i = i + 1

    if found:
        sex = lines[i][:-1]
        l = []
        sexRef = sex[0]
        startRef = 0
        for i in range(1,len(sex)):
            if sex[i] != sexRef:
                l.append((sexRef,startRef+1,i))
                sexRef = sex[i]
                startRef = i
        l.append((sexRef,startRef+1,i+1))
        return l

def get_sexon_id(dictFname):
    '''
    Extrait les identifiants des exons à partir d'un fichier de dictionnaire.
    
    Parameters:
    dictFname (str): Chemin d'accès au fichier de dictionnaire contenant les identifiants des exons.
    
    Returns:
    dict: Un dictionnaire où les clés sont les identifiants d'exons et les valeurs sont les noms correspondants.
    '''
    fdic = open(dictFname)
    lines = fdic.readlines()
    fdic.close()
    d = {}
    for line in lines:
        words = line[:-1].split()
        d[words[1]] = words[0]
    return d



def match_coordinates_and_ids(coordinates_list, id_dict):
    '''Associe les coordonnées des exons avec leurs identifiants.'''
    matched_dict = {}

    for coord in coordinates_list:
        exon_id = id_dict.get(coord[0])
        if exon_id:
            matched_dict[exon_id] = {'start': coord[1], 'end': coord[2]}

    return matched_dict


def calculate_e_value(score, query_length, database_length):
    """
    Calculate the E-value between two sequences.
    
    Args:
    - score: The alignment score between the two sequences.
    - query_length: The length of the query sequence.
    - database_length: The length of the database sequence.
    
    Returns:
    - e_value: The calculated E-value.
    """
    K = 0.041 # K-value, empirically determined for protein sequences
    lambda_value = 0.317 # Lambda value, empirically determined for protein sequences

    # Calculate the raw E-value
    e_value = K * query_length * database_length * math.exp(-lambda_value * score)

    return e_value




def copier_si_inexistant(source, cible):
    '''Copie les éléments du répertoire source vers le 
    répertoire cible s'ils n'existent pas déjà dans ce dernier.'''
    if not os.path.exists(cible):
        os.makedirs(cible)

    for element in os.listdir(source):
        chemin_source = os.path.join(source, element)
        chemin_cible = os.path.join(cible, element)

        if not os.path.exists(chemin_cible):
            shutil.copy(chemin_source, chemin_cible)
            print(f"{element} copié vers {cible}")

def gap(chaine, start, stop):
    '''Ajuste les positions de début et de fin en 
    tenant compte des gaps dans la chaîne.'''
    tiret_indices = []

    for i in range(len(chaine)):
        if chaine[i] == "-":
            tiret_indices.append(i)

    for index in tiret_indices:
        lettres_avant = chaine[:index]
        lettres_apres = chaine[index + 1:] 

        nombre_lettres_avant = len([c for c in lettres_avant if c != "-"])
        nombre_lettres_apres = len([c for c in lettres_apres if c != "-"])

        if nombre_lettres_avant > nombre_lettres_apres:
            stop += 1
        elif nombre_lettres_avant < nombre_lettres_apres:
            start -= 1

    return start, stop

def take_out_dupli(seq):
    '''Supprime les lettres minuscules et les tirets d'une séquence.'''
    new_seq = Seq(''.join(char for char in seq if char.isupper() or char == '-'))
    return new_seq

def calculate_identity_percentage(seq1, seq2):
    '''Calcule le pourcentage d'identité entre deux séquences.'''
    alignments = pairwise2.align.globalxx(seq1, seq2, one_alignment_only=True, score_only=True)
    alignment_length = max(len(seq1), len(seq2))
    identity_percentage = (alignments / alignment_length) * 100
    return identity_percentage


def adding_sequences(sequences_a3m, sequences_msa, exon_start, exon_end, GAP, IDENTITY, exon_id, nouveau_repertoire,nbr_seq):
    # Récupérer la première séquence dans sequences_msa
    first_sequence_msa = sequences_msa[0].seq if sequences_msa else None
    
    for sequence_a3m in sequences_a3m:
        # Sélectionner la séquence entre exon_start et exon_end
        selected_sequence = take_out_dupli(sequence_a3m.seq)[exon_start - 1:exon_end]
        #modifie la séquence si gap
        selected_sequence = gap_inter(first_sequence_msa,selected_sequence)
        # Créer un nouvel objet SeqRecord avec la séquence sélectionnée
        if len(selected_sequence) > 0:
            percentage_gap = (selected_sequence.count("-") / len(selected_sequence)) * 100
            # Ajouter la séquence seulement si le pourcentage de gap est inférieur à 70%
            if percentage_gap <= GAP:
                identity_percentage = calculate_identity_percentage(first_sequence_msa, selected_sequence)
                e_val = calculate_e_value(identity_percentage, len(first_sequence_msa), nbr_seq)
                if e_val < IDENTITY:
                    selected_record = SeqRecord(selected_sequence, id=sequence_a3m.id, description=f"Evalue={e_val}")
                    sequences_msa.append(selected_record)
        else:
            continue
    
    output_msa = nouveau_repertoire + f"msa_s_exon_{exon_id}.fasta"
    SeqIO.write(sequences_msa, output_msa, "fasta")

def verify_exon(exon_id, query_transcript_id, s_exon_table_path):
    # Charger le DataFrame depuis le fichier CSV
    s_exon_table = pd.read_csv(s_exon_table_path)
    
    # Vérifier si l'identifiant d'exon et l'identifiant de transcrit correspondent
    matching_rows = s_exon_table[(s_exon_table['S_exonID'] == exon_id) & (s_exon_table['TranscriptIDCluster'] == query_transcript_id)]
    
    # Vérifier la valeur de la colonne S_exon_Sequence
    if len(matching_rows) == 0:
        # Aucune ligne correspondante trouvée, retourner False
        return False
    else:
        # Au moins une ligne correspondante trouvée, vérifier si S_exon_Sequence est vide
        for index, row in matching_rows.iterrows():
            if pd.isna(row['S_exon_Sequence']):
                # La case est vide, retourner False
                return False
        # Aucune case vide trouvée, retourner True
        return True

def detect_signature(sequence):
    # Créer un dictionnaire avec les positions comme clés et les acides aminés correspondants comme valeurs
    sequence_dict = {}
    for i, aa in enumerate(sequence):
        sequence_dict[i + 1] = aa

    # Définir les positions des acides aminés dans la signature alpha
    positions_alpha = [16, 23, 13, 14, 17, 19, 21]
    # Calculer le nombre d'acides aminés corrects dans la signature alpha
    correct_alpha_count = sum(sequence_dict[pos] == "H" if pos == 16 else
                              sequence_dict[pos] == "R" if pos == 23 else
                              sequence_dict[pos] == "M" if pos == 13 else
                              sequence_dict[pos] == "V" if pos == 14 else
                              sequence_dict[pos] == "K" if pos == 17 else
                              sequence_dict[pos] == "L" if pos == 19 else
                              sequence_dict[pos] == "P" if pos == 21 else False
                              for pos in positions_alpha)
    # Calculer le score de signature alpha
    score_alpha = correct_alpha_count / len(positions_alpha) * 100

    # Définir les positions des acides aminés dans la signature beta
    positions_beta = [16, 23, 11, 15, 18]
    # Calculer le nombre d'acides aminés corrects dans la signature beta
    correct_beta_count = sum(sequence_dict[pos] == "G" if pos == 16 else
                             sequence_dict[pos] == "T" if pos == 23 else
                             sequence_dict[pos] == "G" if pos == 11 else
                             sequence_dict[pos] == "K" if pos == 15 else
                             sequence_dict[pos] == "V" if pos == 18 else False
                             for pos in positions_beta)
    # Calculer le score de signature beta
    score_beta = correct_beta_count / len(positions_beta) * 100

    return score_alpha, score_beta





def adding_sequence_alt(sequences_a3m,sequences_msa,all_msa, positions,exon_start,exon_end, GAP, IDENTITY,exon_id,nouveau_repertoire,nbr_seq,SIGNIFICANT_DIFFERENCE,t,msa_alt_complet):
    exon_id_alt = exon_id + "_alt"
    undecided = []
    
    first_sequence_msa = sequences_msa[0].seq  
    # Initialiser une variable pour stocker la séquence concaténée
    first_sequence_msa_alt = ""

    # Boucler à travers la liste msa_all et concaténer les premières séquences de chaque exon
    for exon_key in all_msa:
        sequences_msa_alt = all_msa[exon_key]
        first_sequence_msa_alt += all_msa[exon_key][0].seq

    ### Build HMM for canonique and alternatif 
    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_msa_file:
        write_msa_to_temp_file(msa_alt_complet, temp_msa_file.name)
        temp_msa_file_path_alt = temp_msa_file.name
    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_msa_file:
        write_msa_to_temp_file(sequences_msa, temp_msa_file.name)
        temp_msa_file_path = temp_msa_file.name

    hmm_file_alt = nouveau_repertoire + f"temp_{exon_id_alt}.hmm"
    hmm_file = nouveau_repertoire + f"temp_{exon_id}.hmm"

    build_hmm(temp_msa_file_path, hmm_file)
    build_hmm(temp_msa_file_path_alt, hmm_file_alt)


    for sequence_a3m in tqdm(sequences_a3m):
        # Sélectionner la séquence entre exon_start et exon_end
        selected_sequence = take_out_dupli(sequence_a3m.seq)[exon_start - 1:exon_end]
        selected_sequence_alt = take_out_dupli(sequence_a3m.seq)[exon_start - 1: exon_start - 1+ 
                                                                 len(first_sequence_msa_alt)]

        selected_sequence = gap_inter(first_sequence_msa,selected_sequence)
        selected_sequence_alt = gap_inter(first_sequence_msa_alt,selected_sequence_alt)
        # Créer un nouvel objet SeqRecord avec la séquence sélectionnée
        if len(selected_sequence) > 0: 
            percentage_gap = (selected_sequence.count("-") / len(selected_sequence)) * 100
            if percentage_gap <= GAP:
                A_score = calculate_identity_percentage(first_sequence_msa, selected_sequence)
                B_score = calculate_identity_percentage(first_sequence_msa_alt, selected_sequence_alt)
                A = calculate_e_value(A_score, len(first_sequence_msa), nbr_seq)
                B = calculate_e_value(B_score, len(first_sequence_msa_alt), nbr_seq)
                Alpha, Beta = detect_signature(selected_sequence)

                if A <= IDENTITY and B <= IDENTITY:
                    A_real_score = score_sequence_with_hmm(hmm_file,str(selected_sequence),exon_id)
                    B_real_score = score_sequence_with_hmm(hmm_file_alt,str(selected_sequence_alt),
                                                            exon_id_alt)
                    pA = math.exp(A_real_score/t)/(math.exp(A_real_score/t)+math.exp(B_real_score/t))
                    pB = math.exp(B_real_score/t)/(math.exp(A_real_score/t)+math.exp(B_real_score/t))
                    os.remove(nouveau_repertoire + f"score_output_{exon_id}.txt")
                    os.remove(nouveau_repertoire + f"score_output_{exon_id_alt}.txt")
                    if pA <= 0.47 or pA >= 0.53:
                        if pA > pB :
                            selected_record = SeqRecord(selected_sequence, id=sequence_a3m.id ,description=f"Evalue={min(A, B)} Alpha={Alpha} Beta={Beta}")
                            sequences_msa.append(selected_record)
                        if pB > pA : 
                            for exon_key in all_msa:
                                begin, end = positions[exon_key]
                                selected_record = SeqRecord(gap_inter(first_sequence_msa_alt, 
                                                                    take_out_dupli(sequence_a3m.seq)[begin - 1:end]), 
                                                                        id=sequence_a3m.id,
                                                                        description=f"Evalue={min(A, B)} Alpha={Alpha} Beta={Beta}")
                                all_msa[exon_key].append(selected_record)
                    else : 
                        undecided_record = SeqRecord(selected_sequence, id=sequence_a3m.id, description=f"Evalue_A={A} Evalue_B={B} Alpha={Alpha} Beta={Beta}")
                        undecided.append(undecided_record)



                #     if Alpha ==  Beta :  
                # # Ajouter la séquence uniquement si l'identité est inférieure ou égale à 60%
                #         if A < B and B - A > SIGNIFICANT_DIFFERENCE:
                #             selected_record = SeqRecord(selected_sequence, id=sequence_a3m.id ,description=f"Evalue={min(A, B)} Alpha={Alpha} Beta={Beta}")
                #             sequences_msa.append(selected_record)
                #         elif B < A and A - B > SIGNIFICANT_DIFFERENCE:
                #             for exon_key in all_msa:
                #                 begin, end = positions[exon_key]
                #                 selected_record = SeqRecord(gap_inter(first_sequence_msa_alt, take_out_dupli(sequence_a3m.seq)[begin - 1:end]), 
                #                                             id=sequence_a3m.id,description=f"Evalue={min(A, B)} Alpha={Alpha} Beta={Beta}")
                #                 all_msa[exon_key].append(selected_record)
                #         else :
                #                 print(sequence_a3m.id)
                #                 undecided_record = SeqRecord(selected_sequence, id=sequence_a3m.id, description=f"Evalue_A={A} Evalue_B={B} Alpha={Alpha} Beta={Beta}")
                #                 undecided.append(undecided_record)
                #     else :
                #         if Alpha > Beta :
                #             selected_record = SeqRecord(selected_sequence, id=sequence_a3m.id ,description=f"Evalue={min(A, B)} Alpha={Alpha} Beta={Beta}")
                #             sequences_msa.append(selected_record)  
                #         else :
                #             for exon_key in all_msa:
                #                 begin, end = positions[exon_key]
                #                 selected_record = SeqRecord(gap_inter(first_sequence_msa_alt, take_out_dupli(sequence_a3m.seq)[begin - 1:end]), 
                #                                                     id=sequence_a3m.id, description=f"Evalue={min(A, B)} Alpha={Alpha} Beta={Beta}")
                #                 all_msa[exon_key].append(selected_record)


    output_msa = nouveau_repertoire + f"msa_s_exon_{exon_id}.fasta"
    SeqIO.write(sequences_msa, output_msa, "fasta")

    output_undecided = inter_path + f'undecided_sequences_{exon_id}.fasta'
    SeqIO.write(undecided, output_undecided, "fasta")
    for exon_key in all_msa : 
        output_msa = nouveau_repertoire + f"msa_s_exon_{exon_key}.fasta"
        SeqIO.write(all_msa[exon_key], output_msa, "fasta")

    # os.remove(nouveau_repertoire + hmm_file)
    # os.remove(nouveau_repertoire + hmm_file_alt)


    

def get_table_asru(ASRU, ID_exons_can):
    df = pd.read_csv(ASRU)
    exon_can_df = pd.DataFrame(columns=range(1))
    exon_alt_df = pd.DataFrame(columns=range(1))

    # Parcourir les lignes de la colonne ASRU
    for _, row in df.iterrows():
        data = row['ASRU'].strip('{}').split(',')  # Séparer par ","
        exon_can = []
        exon_alt = []

        for item in data:
            item = item.strip().strip("'")  # Enlever les apostrophes doubles
            if item in ID_exons_can:
                exon_can.append(item)
            else:
                exon_alt.extend(item.split('$.'))  # Séparer si le symbole "$." est présent

        # Ajouter une nouvelle ligne à la DataFrame pour les exons canoniques
        if exon_can:
            exon_can_df = pd.concat([exon_can_df, pd.DataFrame({0: exon_can})], ignore_index=True)

        # Ajouter une nouvelle ligne à la DataFrame pour les exons alternatifs
        if exon_alt:
            exon_alt_df = pd.concat([exon_alt_df, pd.DataFrame({0: exon_alt})], ignore_index=True)

    return exon_can_df, exon_alt_df



	
def new_borne(exon_start_init, exon_end_init, exon_alt_list):
    all_msa = {}
    positions = {}
    msa_alt_complet = []
    seq_concat = {}  # Dictionnaire pour concaténer les séquences par ID
    first_sequence_msa_alt = ""
    
    for i, exon in enumerate(exon_alt_list):
        sequences_msa = list(SeqIO.parse(msa_directory + f"msa_s_exon_{exon}.fasta", "fasta"))
        all_msa[exon] = sequences_msa
        for seq in sequences_msa:
            if seq.id in seq_concat:
                seq_concat[seq.id] += str(seq.seq)
            else:
                seq_concat[seq.id] = str(seq.seq)
    
    for exon_key in all_msa:
        first_sequence_msa_alt += str(all_msa[exon_key][0].seq)
    
    for i, exon in enumerate(exon_alt_list):
        sequences_msa = list(SeqIO.parse(msa_directory + f"msa_s_exon_{exon}.fasta", "fasta"))
        # Supposons que vous avez une fonction `gap` définie ailleurs
        exon_start, exon_stop = gap(first_sequence_msa_alt, exon_start_init, exon_end_init)
        
        if i == 0:
            start = exon_start
            stop = exon_start + len(sequences_msa[0]) - 1
            stop_prev = stop
        else:
            start = stop_prev + 1
            stop = start + len(sequences_msa[0]) - 1
            stop_prev = stop
        
        positions[exon] = start, stop
    
    # Convertir seq_concat en SeqRecord et l'ajouter à msa_alt_complet
    for id, seq in seq_concat.items():
        msa_alt_complet.append(SeqRecord(Seq(seq), id=id, description=""))
    
    return all_msa, positions, msa_alt_complet




def process_transcript(gene_name, GAP, IDENTITY, SIGNIFICANT_DIFFERENCE,GENE, msa_directory, path_table_path, pir_file_path, dictFname, nouveau_repertoire, ASRU, 
                       transcrit_file,query_transcrit_id,s_exon_table_path,t):
    transcrit_file = pd.read_csv(GENE +'inter/a3m_to_PIR.csv')
    transcript_ids_list = transcrit_file['Transcript IDs'].tolist()
    gene_ids_list = transcrit_file['Gene IDs'].tolist()

    
    index = transcript_ids_list.index(query_transcrit_id)
    query_gene_id = gene_ids_list[index]
    path_table = pd.read_csv(path_table_path)

    if not os.path.exists(nouveau_repertoire):
        os.makedirs(nouveau_repertoire)

    chemin_fichier = os.path.join(nouveau_repertoire, "output.txt")



    coord = get_sexon_coord(query_gene_id, query_transcrit_id, pir_file_path)
    ids = get_sexon_id(dictFname)
    exon_coordinates_transcript = match_coordinates_and_ids(coord, ids)
    csv_file_path = f"DATA/{gene_name}/inter/exon_coordinates_transcript.csv"
    exon_coordinates_df = pd.DataFrame(list(exon_coordinates_transcript.items()), columns=['Key', 'Value'])
    exon_coordinates_df.to_csv(csv_file_path, index=False)
    print(f"Les données ont été enregistrées avec succès dans {csv_file_path}")

    query_exon_paths = path_table.loc[(path_table['GeneID'] == query_gene_id) & (path_table['TranscriptIDCluster'] == query_transcrit_id), 'Path'].tolist()
    exons_similaire = {}
    identifiants = []
    
    for identifiant in query_exon_paths:
        matches = re.findall(r'[^/]+', identifiant)
        identifiants.extend(matches)
    ID_exons_can = set(identifiants)
    asru_table_can, asru_table_alt = get_table_asru(ASRU, ID_exons_can)
    list_empty = []
    asru = []
    for a3m_fichier in glob.glob(GENE + "*.a3m"):
        try:
            sequences_a3m = list(SeqIO.parse(a3m_fichier, "fasta"))[1:]
            nbr_seq = len(sequences_a3m)

            for fichier in tqdm(glob.glob(msa_directory + "msa_s_exon_*.fasta")):
                match = re.search(r"exon_(.*?)\.fasta", fichier)
                if match:
                    exon_id = match.group(1)
                    print("exon : ", exon_id)
                    not_empty = verify_exon(exon_id, query_transcrit_id, s_exon_table_path)
                    if exon_id in ID_exons_can and not not_empty:
                        list_empty.append(exon_id)
                    if exon_id in ID_exons_can and not_empty:
                        exon_start_init = exon_coordinates_transcript[exon_id]['start']
                        exon_end_init = exon_coordinates_transcript[exon_id]['end']
                        sequences_msa = list(SeqIO.parse(fichier, "fasta"))
                        (exon_start, exon_end) = gap(sequences_msa[0].seq, exon_start_init, exon_end_init)

                        if exon_id in asru_table_can.values:

                            print("s-exons similaires detectés",exon_id)
                            index = asru_table_can.columns[asru_table_can.isin([exon_id]).any()].tolist()
                            colonne_index = asru_table_can.columns[asru_table_can.isin([exon_id]).any()]
                            # Convertir le nom de la colonne en index
                            index = asru_table_can.columns.get_loc(colonne_index[0])
                            # Sélectionner la colonne correspondante dans asru_table_alt
                            exon_alt_list = asru_table_alt.iloc[:, index].tolist()
                            L  = [exon_id,exon_alt_list]
                            exons_similaire[exon_id] = exon_alt_list
                            asru.append(L)
                            all_msa, positions,msa_alt_complet = new_borne(exon_start_init, exon_end_init, exon_alt_list)
                            adding_sequence_alt(sequences_a3m,
                                                sequences_msa,
                                                all_msa, positions,
                                                exon_start,
                                                exon_end, GAP, 
                                                IDENTITY,
                                                exon_id,nouveau_repertoire,nbr_seq,SIGNIFICANT_DIFFERENCE,t,msa_alt_complet)
                            


                        else:
                            adding_sequences(sequences_a3m, sequences_msa, 
                                                exon_start, exon_end, GAP,
                                                IDENTITY, exon_id, 
                                                nouveau_repertoire,nbr_seq)

        except ValueError as e:
            print(f"Erreur lors de la lecture du fichier {a3m_fichier}: {e}")

    copier_si_inexistant(msa_directory, nouveau_repertoire)
    with open(chemin_fichier, 'w') as fichier:
        fichier.write(f"gene_name: {gene_name}\n")
        fichier.write(f"query_transcrit_id: {query_transcrit_id}\n")
        fichier.write(f"GAP: {GAP}\n")
        fichier.write(f"IDENTITY: {IDENTITY}\n")
        fichier.write(f"Canonique s-exons empty : {list_empty}\n")
        fichier.write(f"Exons similaires : {asru}\n")
    return query_exon_paths, exons_similaire







def s_exon_augmentation(repertoire, s_exon_table_path):
    output_csv_path = os.path.splitext(s_exon_table_path)[0] + "_a3m.csv"  # Créer le nom de fichier de sortie
    
    with open(s_exon_table_path, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        fieldnames = reader.fieldnames
        
        with open(output_csv_path, 'w', newline='') as output_csvfile:
            writer = csv.DictWriter(output_csvfile, fieldnames=fieldnames)
            writer.writeheader()
            
            for row in reader:
                writer.writerow(row)
    
    with open(output_csv_path, 'a', newline='') as output_csvfile:
        writer = csv.DictWriter(output_csvfile, fieldnames=fieldnames)
        
        for fichier in os.listdir(repertoire):
            if fichier.endswith(".fasta") and fichier.startswith("msa_s_exon_"):
                chemin_fichier = os.path.join(repertoire, fichier)
                identifiant_exon = fichier.split("_exon_")[1].rstrip(".fasta")
                for record in SeqIO.parse(chemin_fichier, "fasta"):
                    nom_transcrit = record.id.split()[0][0:]  # Obtenez l'identifiant du transcrit après ">"
                    if not nom_transcrit.startswith("ENS"):
                        writer.writerow({'TranscriptIDCluster': nom_transcrit,'GeneID': nom_transcrit, 'S_exonID': identifiant_exon})


                        


def rank(table_path, path, exon_similaire):
    # Si path est une liste avec un seul élément, le transformer en chaîne de caractères
    if isinstance(path, list) and len(path) == 1:
        path = path[0]
    
    # Diviser la chaîne de caractères path en une liste d'identifiants d'exon
    path_list = path.split('/')
    
    # Modifier le chemin pour ajouter les exons similaires à la suite de l'exon principal
    for exon_principal, exons_similaires in exon_similaire.items():
        if exon_principal in path_list:
            index_exon_principal = path_list.index(exon_principal)
            path_list = path_list[:index_exon_principal+1] + exons_similaires + path_list[index_exon_principal+1:]
    
    # Lire la table CSV et stocker les données dans une liste de dictionnaires
    with open(table_path, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        table = list(reader)
    
    # Parcourir chaque transcrit ID dans la colonne TranscriptIDCluster
    for transcript_id in set(row['TranscriptIDCluster'] for row in table if not row['TranscriptIDCluster'].startswith("EN")):
        # Initialiser un dictionnaire pour stocker les ordres des exons
        exon_order = {}
        # Initialiser le compteur d'exon
        exon_count = 1
        # Parcourir le chemin et attribuer un numéro unique à chaque exon rencontré
        for exon_id in path_list:
            # Vérifier si l'exon est présent pour le transcrit ID
            if any(row['S_exonID'] == exon_id and row['TranscriptIDCluster'] == transcript_id for row in table):
                # Ajouter l'exon au dictionnaire avec son numéro d'ordre s'il n'est pas déjà présent
                if exon_id not in exon_order:
                    exon_order[exon_id] = exon_count
                    exon_count += 1
        
        # Parcourir chaque ligne de la table pour mettre à jour ExonRank
        for row in table:
            if row['TranscriptIDCluster'] == transcript_id:
                # Récupérer les identifiants d'exon pour ce transcrit ID dans le bon ordre
                exon_ids_ordered = [exon_id for exon_id in path_list if exon_id in row['S_exonID'].split(',')]
                # Créer une liste ordonnée des rangs d'exon
                exon_ranks = [str(exon_order.get(exon_id, '')) for exon_id in exon_ids_ordered]
                # Ajouter les rangs dans la colonne ExonRank
                row['ExonRank'] = ','.join(exon_ranks)

    with open(table_path, 'w', newline='') as output_csvfile:
        writer = csv.DictWriter(output_csvfile, fieldnames=reader.fieldnames)
        writer.writeheader()
        writer.writerows(table)

if __name__ == "__main__":

    if len(sys.argv) != 3:
        print("Usage: python Alignement_ESG.py <gene_name> <transcrit_id>, for the transcrit id please check DATA/gene_name/inter/a3m_to_PIR.csv")
        sys.exit(1)         
    gene_name = sys.argv[1]
    query_transcrit_id = sys.argv[2]

    GAP = 70
    IDENTITY = 1e-4 
    SIGNIFICANT_DIFFERENCE = 1e-5
    t= 100

    GENE = "DATA/" + gene_name + "/"
    msa_directory = GENE + "thoraxe/msa/"
    path_table_path = GENE + "thoraxe/path_table.csv"
    pir_file_path = GENE + 'thoraxe/phylosofs/transcripts.pir'
    dictFname = GENE + '/thoraxe/phylosofs/s_exons.tsv'
    nouveau_repertoire = GENE + "New_alignement/"
    ASRU = GENE + f"antoine_data/{gene_name}_ASRUs_table.csv"
    s_exon_table_path = GENE + "thoraxe/s_exon_table.csv"
    inter_path = GENE + "inter/"
    ases_path = GENE + "thoraxe/ases_table.csv"

    transcrit_file = pd.read_csv(inter_path + 'a3m_to_PIR.csv')

    for a3m_fichier in glob.glob(GENE +"other_data/" + "*.a3m"):
        traiter_fichier_a3m(a3m_fichier)

    exon_path, exon_similaire = process_transcript(gene_name, GAP, IDENTITY,SIGNIFICANT_DIFFERENCE, GENE, msa_directory, path_table_path, pir_file_path, 
                       dictFname, nouveau_repertoire, ASRU, transcrit_file, query_transcrit_id,s_exon_table_path,t)
        

    s_exon_augmentation(nouveau_repertoire,s_exon_table_path)

    new_s_exon_table_path = GENE + "thoraxe/s_exon_table_a3m.csv"
    rank(new_s_exon_table_path,exon_path,exon_similaire)



