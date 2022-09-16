#! /usr/bin/env python
#usage: script.py
import textwrap
import sys
import re

#@profile
def main(argv):
	destfile = sys.argv[2]
	#save contents of read_counts_and_mismatches file as dict per observed species
	#save observed genuses
	#species_counts: {species: [[marker, readcount, correct_bases, total_bases, seqlen, coverage, pid, busco]]}
	species_counts = {}

	counter = 0
	countfile = open(sys.argv[1])
	countfile.readline()
	seq_counts = {}
	species_counts = {}

	for line in countfile:
		counter += 1
		line = line.strip('\n')
		seq = line.split('\t')[1]
		count = int(line.split('\t')[2])
		correct_bases = int(line.split('\t')[3])
		incorrect_bases = int(line.split('\t')[4])
		total_bases = int(line.split('\t')[5])
		subjlen = int(line.split('\t')[6])
		coverage = float(line.split('\t')[7])
		pid = float(line.split('\t')[8])
		org = line.split('\t')[0]

		#save info per sequence in seq_counts dict
		#seq_counts[seq] = [count, correct_bases, total_bases, subjlen, coverage, pid, busco]

		if species not in species_counts:
			species_counts[species] = []
			#find the genus if not a spcollapsed gene

		species_counts[species].append([seq, 
									count, 
									correct_bases, 
									total_bases, 
									subjlen, 
									coverage,
									pid])

	if counter == 0:
		message = "Empty read count file. Likely no aligned reads in sample."
		#print(message)
		#still have to write stuff
		f = open(sys.argv[2], 'w')
		f.write(message + '\n')
		f.close()

		sys.exit()
	countfile.close()

	#done parsing read_counts_and_mismatches file

	#calculate stats for each observed species
	taxon_coverage = {}
	
	#taxon_coverage[taxon] = [observed_markers, 
	#readcounts, 
	#total_bases, 
	#percentage_markers, 
	#marker_coverage, 
	#percent_id, 
	#buscos]

	seen_species = []

	for tax in species_counts:
		mc = len(species_counts[tax])
		counts = 0
		bases = 0
		correct = 0
		total_bases = 0
		subj_len = 0
		buscos = []
		for i in range(0, len(species_counts[tax])):

			busco = species_counts[tax][i][-1]
			if len(busco) > 1:
				buscos.append(busco)

			counts += species_counts[tax][i][1]
			bases += species_counts[tax][i][3]
			correct += species_counts[tax][i][2]
			total_bases += species_counts[tax][i][3]
			subj_len += species_counts[tax][i][4]

		percent_identity = round((correct / total_bases) * 100, 2)
		overall_coverage = round((total_bases / subj_len ) * 100, 2)
		total_markers = len(species_seqs[tax])
		marker_percentage = round( mc / total_markers * 100, 2)
		name = [ncbi.get_species_translator([tax])[e] for e in ncbi.get_species_translator([tax])][0]

		if tax not in seen_speciess:
			seen_speciess.append(tax)

		taxon_coverage[tax] = [mc, 
								counts, 
								total_bases, 
								marker_percentage, 
								overall_coverage, 
								percent_identity, 
								buscos]
	#create tree structure for all observed speciess

	tree = ncbi.get_topology(seen_speciess)
	tree_root = tree.get_tree_root().name
	lineage = ncbi.get_lineage(tree_root)
	tree_speciess = seen_speciess + lineage
	full_tree = ncbi.get_topology(tree_speciess, intermediate_nodes=True)
	full_species_lineage = [node.name for node in full_tree.traverse()]

	#full_seq_speciess: {species: [[specific buscos], specific count, specific + inherited count]}
	full_seq_speciess = {}
	for line in open(files.inherited_markers):
		line = line.strip('\n')
		species = line.split('\t')[0]
		if species in full_species_lineage:
			buscos = []
			for seq in line.split('\t')[1].split(','):
				if len(re.findall('-\d*at\d*-', seq)) > 0:
					busco = re.findall('-\d*at\d*-',seq)[0].strip('-')
					if busco not in buscos:
						buscos.append(busco)

			specific_count = len(line.split('\t')[1].split(','))
			sp_and_inherited_count = len(line.split('\t')[2].split(','))

			full_seq_speciess[species] = [buscos, specific_count, sp_and_inherited_count]

	#determine primary and secondary hits
	#if MRCA is at the level of genus, consider whether one should be primary or secondary by looking at buscos
	primary = {}
	secondary = {}
	for g in genuses:
		if len(genuses[g]) > 1: #multiple species in same genus
			speciess = genuses[g]
			reads = [taxon_coverage[species][1] for species in speciess]
			bases = [taxon_coverage[species][2] for species in speciess]

			#if one has more reads and more bases than all others, it is primary, others are secondary
			maxreads = max(reads)
			maxbases = max(bases)
			pspeciess = []
			if (reads.count(maxreads) == 1 and bases.count(maxbases) == 1)\
			 and (reads.index(maxreads) == bases.index(maxbases)): #no ties, same ID
			 	maxtax = speciess[reads.index(maxreads)]
			 	primary[maxtax] = taxon_coverage[maxtax][0:5]
			 	pspeciess.append(maxtax)
				#pspeciess.append(speciess[reads.index(maxreads)])
				#primary[pspecies] = taxon_coverage[pspecies][0:5]
				#p_buscos = full_seq_speciess[pspecies][0]
			else:
				for t in speciess: 
					if taxon_coverage[t][1] == maxreads or taxon_coverage[t][2] == maxbases:
						pspeciess.append(t)
						primary[t] = taxon_coverage[t][0:5]

			unsorted_aspeciess = [t for t in speciess if t not in pspeciess]
			aspeciess = sorted(unsorted_aspeciess, key = lambda x: taxon_coverage[x][1], reverse = True)
			for aspecies in aspeciess:
				is_secondary = False
				for pspecies in primary:
					p_buscos = full_seq_speciess[pspecies][0]
					a_buscos = taxon_coverage[aspecies][-1]
					a_remain = [b for b in a_buscos if b in p_buscos]
					if len(a_remain) > 0:
						a_above = []
						for b in a_remain:
							#it may not be a hit for the other one! check first
							#check that the pid for this hit is lower
							apid = [seq[6] for seq in species_counts[aspecies] if seq[7] == b]
							ppid = [seq[6] for seq in species_counts[pspecies] if seq[7] == b]
							if len(ppid) > 0 and apid[0] >= ppid[0]:
								a_above.append(b)
							elif len(ppid) == 0:
								a_above.append(b)
						#if a_buscos is fewer than 3, all must be correct
						if len(a_buscos) < 3:
							if len(a_above) < len(a_buscos):
								is_secondary = True
						else:
							if len(a_above) <= len(a_buscos)/3:
								is_secondary = True
					else:
						is_secondary = True
				if is_secondary:
					secondary[aspecies] = taxon_coverage[aspecies][0:5] + [pspecies]
				else:
					primary[aspecies] = taxon_coverage[aspecies][0:5]

		else: #primary
			species = genuses[g][0]
			primary[species] = taxon_coverage[species][0:5]

	#add anything else
	for t in above_species:
		primary[t] = taxon_coverage[t][0:5]

	#write full table
	marker_sorted = sorted(taxon_coverage.keys(), reverse = True, key = lambda x: taxon_coverage[x][3])

	dest = open(files.alltab, 'w')
	dest.write("Name\tTaxid\tObserved_markers\tRead_counts\tPercent_observed_markers\tTotal_marker_coverage\tPercent_identity\n")
	for tax in marker_sorted:
		rank = [ncbi.get_rank([tax])[e] for e in ncbi.get_rank([tax])][0]
		name = [ncbi.get_species_translator([tax])[e] for e in ncbi.get_species_translator([tax])][0]
		mc = taxon_coverage[tax][0]
		counts = taxon_coverage[tax][1]
		marker_percentage = taxon_coverage[tax][3]
		overall_coverage = taxon_coverage[tax][4]
		percent_identity = taxon_coverage[tax][5]
		dest.write(name + '\t'
			+ str(tax) + '\t'
			+ str(mc) + '\t' 
			+ str(counts) + '\t' 
			+ str(marker_percentage) + '%\t'
			+ str(overall_coverage) + '%\t'
			+ str(percent_identity) + '%\n')
	dest.close()

	dest = open(files.primarytab, 'w')
	dest.write("Name\tTaxid\tObserved_markers\tRead_counts\tPercent_observed_markers\tTotal_marker_coverage\tPercent_identity\n")
	#TODO: implement filters

	primary_sorted = sorted(primary.keys(), reverse = True, key = lambda x: primary[x][3])
	#secondary_sorted = sorted(secondary.keys(), reverse=True, key=lambda x: secondary[x][3])

	filter_passing_speciess = []

	for tax in primary_sorted:
		rank = [ncbi.get_rank([tax])[e] for e in ncbi.get_rank([tax])][0]
		name = [ncbi.get_species_translator([tax])[e] for e in ncbi.get_species_translator([tax])][0]
		mc = taxon_coverage[tax][0]
		counts = taxon_coverage[tax][1]
		marker_percentage = taxon_coverage[tax][3]
		overall_coverage = taxon_coverage[tax][4]
		percent_identity = taxon_coverage[tax][5]
		#filter
		if int(mc) >= 2 and int(counts) >= 4:
			filter_passing_speciess.append(tax)
			dest.write(name + '\t'
				+ str(tax) + '\t'
				+ str(mc) + '\t' 
				+ str(counts) + '\t' 
				+ str(marker_percentage) + '%\t'
				+ str(overall_coverage) + '%\t'
				+ str(percent_identity) + '%\n')
	dest.close()

	#close if no filter passing speciess
	if len(filter_passing_speciess) == 0:
		message = "No taxa passing filter requirements."
		#print(message)

		#still have to write stuff
		f = open(files.primarytab, 'w')
		f.write(message + '\n')
		f.close()
		f = open(files.primarytax, 'w')
		f.write(message + '\n')
		f.close()
		sys.exit()

	#create NCBI taxon tree of observed taxa + extend to cellular_org
	tree = ncbi.get_topology(filter_passing_speciess)
	tree_root = tree.get_tree_root().name
	lineage = ncbi.get_lineage(tree_root)
	primary_tree_speciess = [int(e) for e in filter_passing_speciess] + lineage
	primary_tree = ncbi.get_topology(primary_tree_speciess, intermediate_nodes=True)
	#write the tree	structure to file

	orphan_children = []

	#find counts of seqs for internal nodes
	for node in full_tree.traverse():
		if node.is_leaf() == False:
			if node.name not in species_counts:
				species_counts[node.name] = []
			for desc in node.iter_descendants():
				if desc.name in species_counts:
					for seq in species_counts[desc.name]:
						if seq not in species_counts[node.name]:
							species_counts[node.name].append(seq)
		else:
			if node.name not in species_counts:
				orphan_children.append(node.name)


	#create new tree of filter passing hits


	level_counts = []
	currspaces = 0
	currparent = ''
	seen_parents = {}

	dest = open(files.primarytax, 'w')
	dest.write("Markers_Obs\tTotal_Markers\tPercent_Makers_Obs\tPercent_ID\tMarker_read_count\tRank\tName\n")
	for node in primary_tree.traverse("preorder"):

		if node.name not in orphan_children and node.name in full_seq_speciess:
			rank = [ncbi.get_rank([node.name])[e] for e in ncbi.get_rank([node.name])][0]
			name = [ncbi.get_species_translator([node.name])[e] for e in ncbi.get_species_translator([node.name])][0]
			if node.is_root():
				currspaces = 0
			else:
				if currparent == '':
					currparent = node.up.name
					currspaces += 4
				else:
					if currparent != node.up.name:
						currparent = node.up.name
						if currparent in seen_parents:
							currspaces = seen_parents[currparent]
						else:
							currspaces += 4
							seen_parents[currparent] = currspaces
			if node.name in taxon_coverage:
				pid = str(taxon_coverage[node.name][5]) + '%'
			else:
				pid = "NA"
			#total_buscos
			buscos = len(species_counts[str(node.name)])
			seqs = sum([b[1] for b in species_counts[node.name]])
			total_buscos = full_seq_speciess[node.name][2]
			percent = round((buscos/total_buscos)*100,2)
			dest.write(str(buscos) + '\t' 
				+ str(total_buscos) + "\t" 
				+ str(percent)  + '%\t' 
				+ str(pid) + '\t'
				+ str(seqs) + '\t' 
				+ rank + '\t' 
				+ ' ' * currspaces + name + '\n')
	dest.close()

if __name__ == "__main__":
  main(sys.argv)
