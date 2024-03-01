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

def gap_inter(A, B):
    index_tiret = [i+1 for i in range(len(A)-2) if A[i].isalpha() and A[i+1] == '-' and A[i+2].isalpha()]  # Trouver les indices où le motif est trouvé
    for pos in reversed(index_tiret):  # Parcourir en sens inverse pour insérer les tirets correctement
        if pos < len(B):
            B = B[:pos] + '-' + B[pos:]
    return B



# get the exon coordinates from the PIR sequence file
def get_sexon_coord(gid,tid,seqFname):

	# open, read and close the input file
	fseq = open(seqFname)
	lines = fseq.readlines()
	fseq.close()
	# go through the file until finding the tid
	i = 0
	found = False
	while (i < len(lines)) and (not found):
		found = lines[i].startswith('>P1;'+gid+' '+tid)
		i = i + 1
	# if the tid was found
	if found:
		sex = lines[i][:-1]
		# should be a list since they are ordered!
		l = []
		sexRef = sex[0]
		startRef = 0
		# go through the whole sequence
		for i in range(1,len(sex)):
			# if the sexon has just changed
			if sex[i] != sexRef:
				# append the previous sexon to the list
				# +1 to the start to shift
				# nothing to the end because we want the one before last
				l.append((sexRef,startRef+1,i))
				sexRef = sex[i]
				startRef = i
		# don't forget the last one!
		# id only one amino acid, startRef+1 should be the last position
		# and i should be equal to len(sex)
		l.append((sexRef,startRef+1,i+1))
		return l

# get the sexon id from the dictionary file
def get_sexon_id(dictFname):

	# open, read and close the input file
	fdic = open(dictFname)
	lines = fdic.readlines()
	fdic.close()
	d = {}
	for line in lines:
		words = line[:-1].split()
		# key: symbol, value: id
		d[words[1]] = words[0]
	return d

# Match the two last fonctions 
def match_coordinates_and_ids(coordinates_list, id_dict):
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
    # Vérifie si le répertoire cible existe, sinon le crée
    if not os.path.exists(cible):
        os.makedirs(cible)

    # Parcours de tous les éléments du répertoire source
    for element in os.listdir(source):
        chemin_source = os.path.join(source, element)
        chemin_cible = os.path.join(cible, element)

        # Vérifie si l'élément n'existe pas déjà dans le répertoire cible
        if not os.path.exists(chemin_cible):
            # Copie l'élément vers le répertoire cible
            shutil.copy(chemin_source, chemin_cible)
            print(f"{element} copié vers {cible}")

def gap(chaine,start,stop):
    tiret_indices = []

    # Parcourir chaque caractère et trouver les indices des tirets
    for i in range(len(chaine)):
        if chaine[i] == "-":
            tiret_indices.append(i)

    for index in tiret_indices:
        lettres_avant = chaine[:index]


        lettres_apres = chaine[index + 1:] 

        # Compter le nombre de lettres dans lettres_avant (qui ne sont pas "-")
        nombre_lettres_avant = len([c for c in lettres_avant if c != "-"])

        # Compter le nombre de lettres dans lettres_apres (qui ne sont pas "-")
        nombre_lettres_apres = len([c for c in lettres_apres if c != "-"])

        # Comparaison et ajustement de start et stop
        if nombre_lettres_avant > nombre_lettres_apres:
            stop += 1
        elif nombre_lettres_avant < nombre_lettres_apres:
            start -= 1

    return start, stop

def take_out_dupli(seq):
    # Enlever toutes les lettres minuscules de la séquence
    new_seq = Seq(''.join(char for char in seq if char.isupper() or char == '-'))
    return new_seq




def calculate_identity_percentage(seq1, seq2):
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
                    selected_record = SeqRecord(selected_sequence, id=sequence_a3m.id)
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


def adding_sequence_alt(sequences_a3m,sequences_msa,all_msa, positions,exon_start,exon_end, GAP, IDENTITY,exon_id,nouveau_repertoire,nbr_seq,SIGNIFICANT_DIFFERENCE):
    
    first_sequence_msa = sequences_msa[0].seq  
    # Initialiser une variable pour stocker la séquence concaténée
    first_sequence_msa_alt = ""
    
    # Boucler à travers la liste msa_all et concaténer les premières séquences de chaque exon
    for exon_key in all_msa:
        first_sequence_msa_alt += all_msa[exon_key][0].seq
        
        
    for sequence_a3m in sequences_a3m:
        # Sélectionner la séquence entre exon_start et exon_end
        selected_sequence = take_out_dupli(sequence_a3m.seq)[exon_start - 1:exon_end]
        selected_sequence_alt = take_out_dupli(sequence_a3m.seq)[exon_start - 1: exon_start - 1+ len(first_sequence_msa_alt)]

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
                # Ajouter la séquence uniquement si l'identité est inférieure ou égale à 60%
                if A <= IDENTITY and B <= IDENTITY:
                    if A < B and B - A > SIGNIFICANT_DIFFERENCE:
                        selected_record = SeqRecord(selected_sequence, id=sequence_a3m.id)
                        sequences_msa.append(selected_record)
                    elif B < A and A - B > SIGNIFICANT_DIFFERENCE:
                        for exon_key in all_msa : 
                            begin, end = positions[exon_key]
                            selected_record = SeqRecord(gap_inter(first_sequence_msa_alt,take_out_dupli(sequence_a3m.seq)[begin - 1:end]), id=sequence_a3m.id)
                            all_msa[exon_key].append(selected_record)
                    else:
                        continue

    output_msa = nouveau_repertoire + f"msa_s_exon_{exon_id}.fasta"
    SeqIO.write(sequences_msa, output_msa, "fasta")
    for exon_key in all_msa : 
        output_msa = nouveau_repertoire + f"msa_s_exon_{exon_key}.fasta"
        SeqIO.write(all_msa[exon_key], output_msa, "fasta")




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
    all_msa = {}  # Dictionnaire pour stocker les séquences MSA de chaque exon
    positions = {}  # Dictionnaire pour stocker les positions de chaque exon
    first_sequence_msa_alt = ""
    for i, exon in enumerate(exon_alt_list):
        sequences_msa = list(SeqIO.parse(msa_directory + f"msa_s_exon_{exon}.fasta", "fasta"))
        all_msa[exon] = sequences_msa
    
    for exon_key in all_msa:
        first_sequence_msa_alt += all_msa[exon_key][0].seq
    
    for i, exon in enumerate(exon_alt_list):
        sequences_msa = list(SeqIO.parse(msa_directory + f"msa_s_exon_{exon}.fasta", "fasta"))
        exon_start,exon_stop = gap(first_sequence_msa_alt,exon_start_init,exon_end_init)
        
        # Calculer le stop pour l'exon précédent
        if i == 0:
            start = exon_start 
            stop = exon_start + len(sequences_msa[0]) -1
            stop_prev = stop
        else:
            start = stop_prev + 1
            stop = start + len(sequences_msa[0]) - 1
            stop_prev = stop
        
        # Stocker les séquences MSA de cet exon
        
        positions[exon] = start, stop 
    
    return all_msa, positions




def process_transcript(gene_name, GAP, IDENTITY, SIGNIFICANT_DIFFERENCE,GENE, msa_directory, path_table_path, pir_file_path, dictFname, nouveau_repertoire, ASRU, 
                       transcrit_file,query_transcrit_id,s_exon_table_path):
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
    csv_file_path = f"{gene_name}/inter/exon_coordinates_transcript.csv"
    exon_coordinates_df = pd.DataFrame(list(exon_coordinates_transcript.items()), columns=['Key', 'Value'])
    exon_coordinates_df.to_csv(csv_file_path, index=False)
    print(f"Les données ont été enregistrées avec succès dans {csv_file_path}")

    query_exon_paths = path_table.loc[(path_table['GeneID'] == query_gene_id) & (path_table['TranscriptIDCluster'] == query_transcrit_id), 'Path'].tolist()

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

            for fichier in glob.glob(msa_directory + "msa_s_exon_*.fasta"):
                match = re.search(r"exon_(.*?)\.fasta", fichier)
                if match:
                    exon_id = match.group(1)
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
                            asru.append(L)
                            all_msa, positions = new_borne(exon_start_init, exon_end_init, exon_alt_list)

                            adding_sequence_alt(sequences_a3m,
                                                sequences_msa,
                                                all_msa, positions,
                                                exon_start,
                                                exon_end, GAP, 
                                                IDENTITY,
                                                exon_id,nouveau_repertoire,nbr_seq,SIGNIFICANT_DIFFERENCE)

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



if __name__ == "__main__":

    if len(sys.argv) != 3:
        print("Usage: python Alignement_ESG.py <gene_name> <transcrit_id>, for the transcrit id please check DATA/gene_name/inter/a3m_to_PIR.csv")
        sys.exit(1)         
    gene_name = sys.argv[1]
    query_transcrit_id = sys.argv[2]

    GAP = 70
    IDENTITY = 1e-4 
    SIGNIFICANT_DIFFERENCE = 1e-5

    GENE = "../DATA/" + gene_name + "/"
    msa_directory = GENE + "thoraxe/msa/"
    path_table_path = GENE + "thoraxe/path_table.csv"
    pir_file_path = GENE + 'thoraxe/phylosofs/transcripts.pir'
    dictFname = GENE + '/thoraxe/phylosofs/s_exons.tsv'
    nouveau_repertoire = GENE + "New_alignement_biga3m/"
    ASRU = GENE + f"antoine_data/{gene_name}_ASRUs_table.csv"
    s_exon_table_path = GENE + "thoraxe/s_exon_table.csv"

    transcrit_file = pd.read_csv(GENE + 'inter/a3m_to_PIR.csv')

    process_transcript(gene_name, GAP, IDENTITY,SIGNIFICANT_DIFFERENCE, GENE, msa_directory, path_table_path, pir_file_path, 
                       dictFname, nouveau_repertoire, ASRU, transcrit_file, query_transcrit_id,s_exon_table_path)
