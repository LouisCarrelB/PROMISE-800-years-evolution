# authors: Elodie Laine, Louis Carrel-Billiard 
# date of creation: April 4 2024
# purpose: compute the species coverage of s-exons and events detected by ThorAxe
# accounting for their evolutionary distances to a reference species (human, by default)
# WARNING: it seems that the lists of CanonicalPathGenes and AlternativePathGenes from the ases output of Thoraxe
# does not exactly match the lists of genes that actually contain the canonical and alternative s-exons, respectively
# I do not know/remember why. Maybe it's a matter of accounting for the junctions?
# It could also be that there is a bug in ThorAxe (!)
# Concretely, in MAPK8, it looks like the worm does not contain any of the subpath in event 1 while in fact
# three worm genes do contain the 3 s-exons 
# to cope with that I implemented 2 flavours of the assessment for the events:
# - 1 based on the ases table
# - 1 based on the s-exon table, looking at the intersection of the sets of species where the s-exons
# in the canonical (resp. alternative) subpaths (without anchors) are present. 
# WARNING: presence of an s-exon in a species does not imply the associated sequence 
# fits well with the other sequences from the other species...



### 
from Bio import Phylo
import argparse
import pandas as pd
from Bio import BiopythonWarning
import warnings
import bisect
import requests
from tqdm import tqdm 
def get_distances(tree_fname, species_ref):
    # parse the tree file and builds the tree
    tree = Phylo.read(tree_fname, "newick")
    
    # get the names of the leaves
    leaves = {term.name for term in tree.get_terminals()}
    
    # Initialiser le dictionnaire pour stocker les distances
    d = {}
    
    # Calculer la distance pour chaque feuille par rapport à species_ref
    for l in leaves:
        # Empiriquement, il semble que nous devrions diviser par 2
        d[l.lower()] = tree.distance(species_ref, l) / 2
    
    return d


def get_species_per_sex(s_exon_fname):

	with open(s_exon_fname) as f:
		lines = f.readlines()

	dsex = {}
	dgene = {}
	for line in lines[1:]:
		wds = line[:-1].split(',')
		species = wds[0]
		geneId = wds[1]
		sexon = wds[5]
		if sexon not in dsex:
			dsex[sexon] = set()
		dsex[sexon].add(species)
		if geneId not in dgene:
			dgene[geneId] = species

	return (dsex,dgene)

def get_species_per_subpath(events_fname, species_d):

	with open(events_fname) as f:
		lines = f.readlines()

	levents = []
	for line in lines[1:]:
		# CanonicalPath: 0, AlternativePath: 1
		wds = line[:-1].split(',')
		CanonicalPath = wds[0].split("/")[1:-1]
		AlternativePath = wds[1].split("/")[1:-1]
		labels = ['can','alt']
		s_exons = (CanonicalPath, AlternativePath)
		d = {}	
		for i in range(2):
			if len(s_exons[i])>0:
				d[labels[i]] = species_d[s_exons[i][0]]
				for s in s_exons[i][1:]:
					d[labels[i]] = d[labels[i]].intersection(species_d[s])
			else:
				d[labels[i]] = set()
		levents.append(d)
	#print(levents[1:3])
	return (levents)

# count the number of species in each predefined evolutionary distance bin
# the output are 2 dictionaries
# one with key: species, value: associated bin (lower than that evol dist) 
# the other one with key: evol dist and value: count
# dist_d is a dictionary with species as keys and evol distances as values
def discretize(dist_d, bins=[10,50,200,400,600,800,1000,1500]):
	counts = {}
	d_bins = {}
	# for each species
	for s in dist_d:
		# get its bin index
		index = bisect.bisect_right(bins, dist_d[s])
		# associate the species with the bin
		d_bins[s] = bins[index]
		# count +1 for this bin
		if bins[index] not in counts:
			counts[bins[index]] = 0
		counts[bins[index]] += 1
	print(dict(sorted(counts.items())))
	return d_bins, counts
 
# species_d is a dictionary with s-exons as keys and list of species as values
def calc_evol_cov_s_exons(species_d, bins_d, counts_tot) :
	sex_d ={}
	s_missing = set()
	# for each s-exon
	for sex in species_d:
		# add the s-exon if it's not already there
		if sex not in sex_d:
			sex_d[sex] = {}
		# for each species where that s-exon is present
		for species in species_d[sex]:
			# if the species is associated with a bin			
			if species in bins_d:
				# count +1 in that bin for the current s-exon
				if bins_d[species] not in sex_d[sex]:
					sex_d[sex][bins_d[species]] = 0
				sex_d[sex][bins_d[species]] += 1.0/counts_tot[bins_d[species]]
			else:
				s_missing.add(species)
	print("species not found for s-exons:", s_missing)
	return sex_d


def calc_evol_cov_events(events_fname, species_gene_d, bins_d, counts_tot):

	with open(events_fname) as f:
		lines = f.readlines()

	s_missing = set()
	levents = []
	for line in lines[1:]:
		# CanonicalPathGenes: 8, AlternativePathGenes: 9
		wds = line[:-1].split(',')
		CanonicalPathGenes = wds[8].split("/")
		AlternativePathGenes = wds[9].split("/")
		labels = ['can','alt']
		genes = (CanonicalPathGenes, AlternativePathGenes)	
		d = {}
		for i in range(2):
			# initialize the dictionary of the can or alt subpath
			d[labels[i]] = {}
			species_visited = []
			# for each species where the can or alt subpath is expressed
			for gene in genes[i]:
				species = species_gene_d[gene]
				if species not in species_visited:
					if species in bins_d:
						if bins_d[species] not in d[labels[i]]:
							d[labels[i]][bins_d[species]] = 0
						d[labels[i]][bins_d[species]] += 1.0/counts_tot[bins_d[species]]
					else:
						s_missing.add(species)
					species_visited.append(species)
		levents.append(d)
	print("species not found for events:", s_missing)
	return levents

def calc_evol_cov_events2(levents, bins_d, counts_tot):

	labels = ['can','alt']
	s_missing = set()
	res = []
	for event in levents:
		d = {}
		for i in range(2):
			# initialize the dictionary of the can or alt subpath
			d[labels[i]] = {}
			# for each species where the can or alt subpath is expressed
			for species in event[labels[i]]:
				if species in bins_d:
					if bins_d[species] not in d[labels[i]]:
						d[labels[i]][bins_d[species]] = 0
					d[labels[i]][bins_d[species]] += 1.0/counts_tot[bins_d[species]]
				else:
					s_missing.add(species)
		res.append(d)
	print("species not found for events:", s_missing)
	return res

def write_results_sex(sex_d, fname, bins=[10,50,200,400,600,800,1000,1500]):
	with open(fname,'w') as f:
		f.write('Sexon_ID,Dist,Fraction\n')
		for sex in sex_d:
			for dist in bins:
				if dist in sex_d[sex]:
					f.write(sex+','+str(dist)+','+ástr(sex_d[sex][dist])+'\n')
				else:
					f.write(sex+','+str(dist)+',0\n')

def write_results_event(event_l, fname, bins=[10,50,200,400,600,800,1000,1500]):
	with open(fname,'w') as f:
		f.write('Event,Subpath,Dist,Fraction\n')
		labels = ['can','alt']
		for i in range(len(event_l)):
			for k in range(2):
				for dist in bins:
					if dist in event_l[i][labels[k]]:
						f.write(str(i+1)+','+labels[k]+','+str(dist)+','+str(event_l[i][labels[k]][dist])+'\n')
					else:
						f.write(str(i+1)+','+labels[k]+','+str(dist)+',0\n')


def standardize_species_names(file_path):
    # Charger le fichier CSV
    df = pd.read_csv(file_path)
    
    # Vérifier si la colonne "Species" existe
    if 'Species' in df.columns:
        # Standardiser les noms des espèces
        df['Species'] = df['Species'].apply(lambda x: x.lower().replace(' ', '_'))
        
        # Enregistrer le fichier modifié
        df.to_csv(file_path, index=False)
        print("Les noms des espèces ont été standardisés.")
    else:
        print("La colonne 'Species' n'existe pas dans le fichier.")



def adjust_index(row):
    if row['Ambigue'] == 'True':
        if row['Index'] == 'Not_in_path':
            return 'Ambiguous'
        elif row['Index'] in ['can', 'alt', 'both']:
            return row['Index'] + '_and_ambi'
    return row['Index']



def capitalize_and_format(species_name):
    if species_name == "":
        return ""  # Retourner une chaîne vide si l'entrée est vide
    if species_name.lower() == "blastocystis hominis":
        formatted_name = species_name.replace(' ', '_').capitalize()
    else:
        formatted_name = species_name[0].upper() + species_name[1:]
    return formatted_name


def format_species_names(input_file_path, output_file_path):
    # Charger le fichier CSV avec la détection de la tabulation comme séparateur
    df = pd.read_csv(input_file_path, sep='\t')
    
    # Transformer les noms des espèces
    df['Species'] = df['Species'].apply(lambda x: ' '.join(word.capitalize() for word in x.split('_')))
    
    # Créer un nouveau DataFrame contenant uniquement la colonne "Species"
    new_df = df[['Species']]
    
    # Enregistrer le nouveau DataFrame dans un fichier CSV sans la colonne "Index"
    new_df.to_csv(output_file_path, index=False)




def get_uniref_id(full_id):
    """Extrait l'identifiant UniProt depuis un identifiant UniRef complet."""
    return full_id.split('_')[1] if '_' in full_id else full_id

def get_organism_name(organism_info):
    names = organism_info.get("names", [])
    for name_entry in names:
        if name_entry.get("type", "") == "scientific":
            return name_entry.get("value", "")
    return None

def get_organism_info(uniprot_id):
    if uniprot_id.startswith("Uni"):
        url = f"https://www.ebi.ac.uk/proteins/api/proteins/{get_uniref_id(uniprot_id)}"
    else:
        url = f"https://www.ebi.ac.uk/proteins/api/proteins/{uniprot_id}"
    response = requests.get(url)
    if response.status_code == 200:
        protein_data = response.json()
        organism_info = protein_data.get("organism", {})
        organism_name = get_organism_name(organism_info)
        return organism_name
    return None

def read_msa_and_get_organisms(filename):
    organism_list = []
    with open(filename, 'r') as file:
        for line in tqdm(file, desc=  "Looking for species in undecided sequences"):
            if line.startswith('>'):
                uniprot_id = line.split()[0][1:]  # Assume ID follows '>'
                organism_name = get_organism_info(uniprot_id)
                if organism_name:
                    organism_list.append(organism_name)
    return organism_list


if __name__ == "__main__":

	def arg_parser():  # Parser implementation
		parser = argparse.ArgumentParser(prog='get_evol_dist_coverage.py',
											epilog='get_evol_dist_coverage.py --gid ENSG00000107643 --species Homo_sapiens',
											formatter_class=argparse.RawDescriptionHelpFormatter)
		parser.add_argument('--gid', help='Ensembl gene identifier, eg: ENSG00000107643')
		parser.add_argument('--species', help='Reference species, eg: Homo_sapiens')
		return parser.parse_args()

	# parse and config arguments
	args = arg_parser()
	if args.species is None:
		args.species = 'Homo_sapiens' 

	tree_fname = "DATA/" + args.gid+'/list_species.nwk.txt'
	s_exon_fname = "DATA/" +args.gid+'/thoraxe/s_exon_table_a3m.csv'
	events_fname = "DATA/" +args.gid+'/thoraxe/ases_table_a3m.csv'
	bins = [10,50,200,400,600,800,1000,1500]
	# Utiliser la fonction
	standardize_species_names(s_exon_fname)
	# # get the evolutionary distances of each species wrt the reference species (human by default)
	# ² = get_distances(tree_fname, args.species)
	# # print(dist_d)
	# max_val = max(dist_d.values())
	# print('Max evolutionary distance:', max_val)
	# bins_d, counts_tot = discretize(dist_d, bins)
	# # get the list of species where each s-exon is present
	species_d, species_gene_d = get_species_per_sex(s_exon_fname)
	
	

	# # s-exon coverage
	# sex_d = calc_evol_cov_s_exons(species_d, bins_d, counts_tot)
	
	# fout = args.gid + '/thoraxe_2/s_exon_evol_cov.csv'
	# write_results_sex(sex_d, fout, bins)

	# # event coverage from the geneIds indicated in ThorAxe output
	# event_l = calc_evol_cov_events(events_fname, species_gene_d, bins_d, counts_tot)
	
	# fout = args.gid + '/thoraxe_2/events_evol_cov.csv'
	# write_results_event(event_l, fout, bins)

		# event coverage by directly looking at the representation of s-exons sets (subpaths)
	event_l = get_species_per_subpath(events_fname, species_d)
	firstevent = event_l[0]
	

	# event_l = calc_evol_cov_events2(event_l, bins_d, counts_tot)
	# print(event_l)
	# fout = args.gid + '/thoraxe_2/events_evol_cov_from_sex_table.csv'
	# write_results_event(event_l, fout, bins)
	# df = pd.read_csv(s_exon_fname)
	# print(df.columns)  # Ceci affichera les noms des colonnes dans la console

	
	
	
############# ADDING FROM LOUIS ASES
	# Créer un dictionnaire pour enregistrer les clés associées à chaque espèce
	## Ambigues :
	filename = "DATA/" +args.gid+"/inter/undecided_sequences_17_0.fasta"
	ambigue = read_msa_and_get_organisms(filename)
	ambigue = [s.lower().replace(' ', '_') for s in ambigue]



	species_keys = {}
	output_filename = "DATA/" +args.gid + '/filtered_species_for_pastml.tsv'
	# Itérer sur chaque paire clé-ensemble du dictionnaire
	for key, species_set in firstevent.items():
		for species in species_set:
			if species in species_keys:
				species_keys[species].add(key)
			else:
				species_keys[species] = {key}

	# Créer une liste pour stocker les tuples (index, espèce)
	data = []

	# Itérer sur le dictionnaire species_keys pour créer les données finales
	for species, keys in species_keys.items():
		index_value = 'both' if len(keys) > 1 else list(keys)[0]
		data.append((index_value, species))

	# Convertir en DataFrame
	df = pd.DataFrame(data, columns=['Index', 'Species'])

	# Présumant que args.gid est une variable contenant le chemin du dossier
	gid_path = args.gid  # Remplacez ceci par la valeur réelle de args.gid
	output_file_path = f"DATA/{gid_path}/filtered_species_for_pastml.csv"

	
	species = set()
	for value_set in species_d.values():
		species.update(value_set)
	existing_species = set(df['Species'])

	# Déterminer les espèces manquantes
	missing_species = species - existing_species
	print(species)
	# Créer les nouvelles lignes pour les espèces manquantes
	new_rows = [{'Index': 'Not_in_path', 'Species': sp} for sp in missing_species]

	# Ajouter les nouvelles lignes au DataFrame existant
	df = pd.concat([df, pd.DataFrame(new_rows)], ignore_index=True)
	



	df['Ambigue'] = df['Species'].apply(lambda x: 'True' if x in ambigue else 'False')
	ambigue_species_not_in_df = set(ambigue) - set(df['Species'])
	new_rows_ambigue = [{'Index': 'Not_in_path', 'Species': sp, 'Ambigue': 'True'} for sp in ambigue_species_not_in_df]
	df = pd.concat([df, pd.DataFrame(new_rows_ambigue)], ignore_index=True)

	df['Index'] = df.apply(adjust_index, axis=1)

	df = df.drop(columns=['Ambigue'])
	print(df)
	



	df['Species'] = df['Species'].apply(capitalize_and_format)
	

	df.to_csv(output_filename, sep='\t', index=False, columns=['Species', 'Index'])

	species_path =  "DATA/" + args.gid + '/Species_list_time_tree.csv'
	format_species_names(output_filename, species_path)


	




