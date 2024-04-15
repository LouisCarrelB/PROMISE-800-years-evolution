# author: Elodie Laine
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

from Bio import Phylo
import argparse
import pandas as pd
from Bio import BiopythonWarning
import warnings
import bisect

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
					f.write(sex+','+str(dist)+','+str(sex_d[sex][dist])+'\n')
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

    tree_fname = args.gid+'/list_species.nwk.txt'
    s_exon_fname = args.gid+'/thoraxe/s_exon_table.csv'
    events_fname = args.gid+'/thoraxe/ases_table.csv'

    bins=[10,50,200,400,600,800,1000,1500]

    # get the evolutionary distances of each species wrt the reference species (human by default)
    dist_d = get_distances(tree_fname,args.species)
    #print(dist_d)
    max_val = max(dist_d.values())
    print('Max evolutionary distance:', max_val)
    bins_d, counts_tot = discretize(dist_d,bins)
    # get the list of species where each s-exon is present
    species_d, species_gene_d = get_species_per_sex(s_exon_fname)

    # s-exon coverage
    sex_d = calc_evol_cov_s_exons(species_d, bins_d, counts_tot)
    fout = args.gid+'/thoraxe/s_exon_evol_cov.csv'
    write_results_sex(sex_d, fout, bins)

    # event coverage from the geneIds indicated in ThorAxe output
    event_l = calc_evol_cov_events(events_fname, species_gene_d, bins_d, counts_tot)
    fout = args.gid+'/thoraxe/events_evol_cov.csv'
    write_results_event(event_l, fout, bins)

    # event coverage by directly looking at the representation of s-exons sets (subpaths)
    event_l = get_species_per_subpath(events_fname, species_d)
    event_l = calc_evol_cov_events2(event_l, bins_d, counts_tot)
    fout = args.gid+'/thoraxe/events_evol_cov_from_sex_table.csv'
    write_results_event(event_l, fout, bins)
    
