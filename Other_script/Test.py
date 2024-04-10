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



